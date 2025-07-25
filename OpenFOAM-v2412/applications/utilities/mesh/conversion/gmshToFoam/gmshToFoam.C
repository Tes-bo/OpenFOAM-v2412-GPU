/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    gmshToFoam

group
    grpMeshConversionUtilities

Description
    Reads .msh file as written by Gmsh.

    Needs surface elements on mesh to be present and aligned with outside faces
    of the mesh. I.e. if the mesh is hexes, the outside faces need to be
    quads.

    Note: There is something seriously wrong with the ordering written in the
    .msh file. Normal operation is to check the ordering and invert prisms
    and hexes if found to be wrong way round.
    Use the -keepOrientation to keep the raw information.

    Note: The code now uses the element (cell,face) physical region id number
    to create cell zones and faces zones (similar to
    fluentMeshWithInternalFaces).

    A use of the cell zone information, is for field initialization with the
    "setFields" utility. see the classes:  topoSetSource, zoneToCell.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "cellModel.H"
#include "repatchPolyTopoChanger.H"
#include "cellSet.H"
#include "faceSet.H"
#include "List.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Element type numbers

static label MSHLINE   = 1;

static label MSHTRI   = 2;
static label MSHQUAD  = 3;
static label MSHTET   = 4;


static label MSHHEX   = 5;
static label MSHPRISM = 6;
static label MSHPYR   = 7;


// Skips till end of section. Returns false if end of file.
bool skipSection(IFstream& inFile)
{
    string line;
    do
    {
        inFile.getLine(line);

        if (!inFile.good())
        {
            return false;
        }
    }
    while (line.size() < 4 || line.substr(0, 4) != "$End");

    return true;
}


// Find face in pp which uses all vertices in meshF (in mesh point labels)
label findFace(const primitivePatch& pp, const labelList& meshF)
{
    const Map<label>& meshPointMap = pp.meshPointMap();

    // meshF[0] in pp labels.
    if (!meshPointMap.found(meshF[0]))
    {
        Warning<< "Not using gmsh face " << meshF
            << " since zero vertex is not on boundary of polyMesh" << endl;
        return -1;
    }

    // Find faces using first point
    const labelList& pFaces = pp.pointFaces()[meshPointMap[meshF[0]]];

    // Go through all these faces and check if there is one which uses all of
    // meshF vertices (in any order ;-)
    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        const face& f = pp[facei];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (meshF.found(f[fp]))
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return facei;
        }
    }

    return -1;
}


// Same but find internal face. Expensive addressing.
label findInternalFace(const primitiveMesh& mesh, const labelList& meshF)
{
    const labelList& pFaces = mesh.pointFaces()[meshF[0]];

    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        const face& f = mesh.faces()[facei];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (meshF.found(f[fp]))
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return facei;
        }
    }
    return -1;
}


// Determine whether cell is inside-out by checking for any wrong-oriented
// face.
bool correctOrientation(const pointField& points, const cellShape& shape)
{
    // Get centre of shape.
    const point cc(shape.centre(points));

    // Get outwards pointing faces.
    faceList faces(shape.faces());

    for (const face& f : faces)
    {
        const vector areaNorm(f.areaNormal(points));

        // Check if vector from any point on face to cc points outwards
        if (((points[f[0]] - cc) & areaNorm) < 0)
        {
            // Incorrectly oriented
            return false;
        }
    }

    return true;
}


void storeCellInZone
(
    const label regPhys,
    const label celli,
    Map<label>& physToZone,

    labelList& zoneToPhys,
    List<DynamicList<label>>& zoneCells
)
{
    const auto zoneFnd = physToZone.cfind(regPhys);

    if (zoneFnd.good())
    {
        // Existing zone for region
        zoneCells[zoneFnd()].append(celli);
    }
    else
    {
        // New region. Allocate zone for it.
        const label zonei = zoneCells.size();
        zoneCells.setSize(zonei+1);
        zoneToPhys.setSize(zonei+1);

        Info<< "Mapping region " << regPhys << " to Foam cellZone "
            << zonei << endl;
        physToZone.insert(regPhys, zonei);

        zoneToPhys[zonei] = regPhys;
        zoneCells[zonei].append(celli);
    }
}


// Reads mesh format
scalar readMeshFormat(IFstream& inFile)
{
    Info<< "Starting to read mesh format at line "
        << inFile.lineNumber()
        << endl;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    scalar version;
    label asciiFlag, nBytes;
    lineStr >> version >> asciiFlag >> nBytes;

    Info<< "Read format version " << version << "  ascii " << asciiFlag << endl;

    if (asciiFlag != 0)
    {
        FatalIOErrorInFunction(inFile)
            << "Can only read ascii msh files."
            << exit(FatalIOError);
    }

    inFile.getLine(line);
    IStringStream tagStr(line);
    word tag(tagStr);

    if (tag != "$EndMeshFormat")
    {
        FatalIOErrorInFunction(inFile)
            << "Did not find $ENDNOD tag on line "
            << inFile.lineNumber() << exit(FatalIOError);
    }
    Info<< endl;

    return version;
}


// Reads points and map for gmsh MSH file <4
void readPointsLegacy(IFstream& inFile, pointField& points, Map<label>& mshToFoam)
{
    Info<< "Starting to read points at line " << inFile.lineNumber() << endl;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    label nVerts;
    lineStr >> nVerts;

    Info<< "Vertices to be read: " << nVerts << endl;

    points.resize(nVerts);

    for (label pointi = 0; pointi < nVerts; pointi++)
    {
        label mshLabel;
        scalar xVal, yVal, zVal;

        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);

        lineStr >> mshLabel >> xVal >> yVal >> zVal;

        point& pt = points[pointi];

        pt.x() = xVal;
        pt.y() = yVal;
        pt.z() = zVal;

        mshToFoam.insert(mshLabel, pointi);
    }

    Info<< "Vertices read:" << mshToFoam.size() << endl;

    inFile.getLine(line);
    IStringStream tagStr(line);
    word tag(tagStr);

    if (tag != "$ENDNOD" && tag != "$EndNodes")
    {
        FatalIOErrorInFunction(inFile)
            << "Did not find $ENDNOD tag on line "
            << inFile.lineNumber() << exit(FatalIOError);
    }
    Info<< endl;
}

// Reads points and map for gmsh MSH file >=4
void readPoints(IFstream& inFile, pointField& points, Map<label>& mshToFoam)
{
    Info<< "Starting to read points at line " << inFile.lineNumber() << endl;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    // Number of "entities": 0, 1, 2, and 3 dimensional geometric structures
    label nEntities, nVerts;
    lineStr >> nEntities >> nVerts;

    Info<< "Vertices to be read: " << nVerts << endl;

    points.resize(nVerts);

    // Index of points, in the order as they appeared, not what gmsh
    // labelled them.
    label pointi = 0;

    for (label entityi = 0; entityi < nEntities; entityi++)
    {
        label entityDim, entityLabel, isParametric, nNodes;
        scalar xVal, yVal, zVal;
        inFile.getLine(line);
        IStringStream lineStr(line); // can IStringStream be removed?

        // Read entity entry, then set up a list for node IDs
        lineStr >> entityDim >> entityLabel >> isParametric >> nNodes;
        List<label> nodeIDs(nNodes);

        // Loop over entity node IDs
        for (label eNode = 0; eNode < nNodes; ++eNode)
        {
            inFile.getLine(line);
            IStringStream lineStr(line);
            lineStr >> nodeIDs[eNode];
        }

        // Loop over entity node values, saving to points[]
        for (label eNode = 0; eNode < nNodes; ++eNode)
        {
            inFile.getLine(line);
            IStringStream lineStr(line);
            lineStr >> xVal >> yVal >> zVal;
            point& pt = points[nodeIDs[eNode]-1];
            pt.x() = xVal;
            pt.y() = yVal;
            pt.z() = zVal;
            mshToFoam.insert(nodeIDs[eNode], pointi++);
        }

    }

    Info<< "Vertices read: " << mshToFoam.size() << endl;

    inFile.getLine(line);
    IStringStream tagStr(line);
    word tag(tagStr);

    if (tag != "$ENDNOD" && tag != "$EndNodes")
    {
        FatalIOErrorInFunction(inFile)
            << "Did not find $ENDNOD tag on line "
            << inFile.lineNumber() << exit(FatalIOError);
    }
    Info<< endl;
}


// Reads physical names
void readPhysNames(IFstream& inFile, Map<word>& physicalNames)
{
    Info<< "Starting to read physical names at line " << inFile.lineNumber()
        << endl;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    label nNames;
    lineStr >> nNames;

    Info<< "Physical names:" << nNames << endl;

    for (label i = 0; i < nNames; i++)
    {
        label regionI;
        string regionName;

        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);
        label nSpaces = lineStr.str().count(' ');

        if (nSpaces == 1)
        {
            lineStr >> regionI >> regionName;

            Info<< "    " << regionI << '\t'
                << word::validate(regionName) << endl;
        }
        else if (nSpaces == 2)
        {
            // >= Gmsh2.4 physical types has tag in front.
            label physType;
            lineStr >> physType >> regionI >> regionName;
            if (physType == 1)
            {
                Info<< "    " << "Line " << regionI << '\t'
                    << word::validate(regionName) << endl;
            }
            else if (physType == 2)
            {
                Info<< "    " << "Surface " << regionI << '\t'
                    << word::validate(regionName) << endl;
            }
            else if (physType == 3)
            {
                Info<< "    " << "Volume " << regionI << '\t'
                    << word::validate(regionName) << endl;
            }
        }
        else
        {
            continue;
        }

        physicalNames.insert(regionI, word::validate(regionName));
    }

    inFile.getLine(line);
    IStringStream tagStr(line);
    word tag(tagStr);

    if (tag != "$EndPhysicalNames")
    {
        FatalIOErrorInFunction(inFile)
            << "Did not find $EndPhysicalNames tag on line "
            << inFile.lineNumber() << exit(FatalIOError);
    }
    Info<< endl;
}

void readCellsLegacy
(
    const scalar versionFormat,
    const bool keepOrientation,
    const pointField& points,
    const Map<label>& mshToFoam,
    IFstream& inFile,
    cellShapeList& cells,

    labelList& patchToPhys,
    List<DynamicList<face>>& patchFaces,

    labelList& zoneToPhys,
    List<DynamicList<label>>& zoneCells
)
{
    //$Elements
    //number-of-elements
    //elm-number elm-type number-of-tags < tag > \u2026 node-number-list


    Info<< "Starting to read cells at line " << inFile.lineNumber() << endl;

    const cellModel& hex = cellModel::ref(cellModel::HEX);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& pyr = cellModel::ref(cellModel::PYR);
    const cellModel& tet = cellModel::ref(cellModel::TET);

    face triPoints(3);
    face quadPoints(4);
    labelList tetPoints(4);
    labelList pyrPoints(5);
    labelList prismPoints(6);
    labelList hexPoints(8);


    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    label nElems;
    lineStr >> nElems;

    Info<< "Cells to be read: " << nElems << endl << endl;


    // Storage for all cells. Too big. Shrink later
    cells.setSize(nElems);

    label celli = 0;
    label nTet = 0;
    label nPyr = 0;
    label nPrism = 0;
    label nHex = 0;


    // From gmsh physical region to Foam patch
    Map<label> physToPatch;

    // From gmsh physical region to Foam cellZone
    Map<label> physToZone;


    for (label elemI = 0; elemI < nElems; elemI++)
    {
        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);

        label elmNumber, elmType, regPhys;
        if (versionFormat >= 2)
        {
            lineStr >> elmNumber >> elmType;

            label nTags;
            lineStr >> nTags;

            if (nTags > 0)
            {
                // Assume the first tag is the physical surface
                lineStr >> regPhys;
                for (label i = 1; i < nTags; i++)
                {
                    label dummy;
                    lineStr >> dummy;
                }
            }
        }
        else
        {
            label regElem, nNodes;
            lineStr >> elmNumber >> elmType >> regPhys >> regElem >> nNodes;
        }

        // regPhys on surface elements is region number.
        if (elmType == MSHLINE)
        {
            label meshPti;
            lineStr >> meshPti >> meshPti;
        }
        else if (elmType == MSHTRI)
        {
            lineStr >> triPoints[0] >> triPoints[1] >> triPoints[2];

            inplaceRenumber(mshToFoam, triPoints);

            const auto regFnd = physToPatch.cfind(regPhys);

            label patchi = -1;
            if (regFnd.good())
            {
                // Existing patch for region
                patchi = regFnd();
            }
            else
            {
                // New region. Allocate patch for it.
                patchi = patchFaces.size();

                patchFaces.setSize(patchi + 1);
                patchToPhys.setSize(patchi + 1);

                Info<< "Mapping region " << regPhys << " to Foam patch "
                    << patchi << endl;
                physToPatch.insert(regPhys, patchi);
                patchToPhys[patchi] = regPhys;
            }

            // Add triangle to correct patchFaces.
            patchFaces[patchi].append(triPoints);
        }
        else if (elmType == MSHQUAD)
        {
            lineStr
                >> quadPoints[0] >> quadPoints[1] >> quadPoints[2]
                >> quadPoints[3];

            inplaceRenumber(mshToFoam, quadPoints);

            const auto regFnd = physToPatch.cfind(regPhys);

            label patchi = -1;
            if (regFnd.good())
            {
                // Existing patch for region
                patchi = regFnd();
            }
            else
            {
                // New region. Allocate patch for it.
                patchi = patchFaces.size();

                patchFaces.setSize(patchi + 1);
                patchToPhys.setSize(patchi + 1);

                Info<< "Mapping region " << regPhys << " to Foam patch "
                    << patchi << endl;
                physToPatch.insert(regPhys, patchi);
                patchToPhys[patchi] = regPhys;
            }

            // Add quad to correct patchFaces.
            patchFaces[patchi].append(quadPoints);
        }
        else if (elmType == MSHTET)
        {
            storeCellInZone
            (
                regPhys,
                celli,
                physToZone,
                zoneToPhys,
                zoneCells
            );

            lineStr
                >> tetPoints[0] >> tetPoints[1] >> tetPoints[2]
                >> tetPoints[3];

            inplaceRenumber(mshToFoam, tetPoints);

            cells[celli++].reset(tet, tetPoints);

            nTet++;
        }
        else if (elmType == MSHPYR)
        {
            storeCellInZone
            (
                regPhys,
                celli,
                physToZone,
                zoneToPhys,
                zoneCells
            );

            lineStr
                >> pyrPoints[0] >> pyrPoints[1] >> pyrPoints[2]
                >> pyrPoints[3] >> pyrPoints[4];

            inplaceRenumber(mshToFoam, pyrPoints);

            cells[celli++].reset(pyr, pyrPoints);

            nPyr++;
        }
        else if (elmType == MSHPRISM)
        {
            storeCellInZone
            (
                regPhys,
                celli,
                physToZone,
                zoneToPhys,
                zoneCells
            );

            lineStr
                >> prismPoints[0] >> prismPoints[1] >> prismPoints[2]
                >> prismPoints[3] >> prismPoints[4] >> prismPoints[5];

            inplaceRenumber(mshToFoam, prismPoints);

            cells[celli].reset(prism, prismPoints);

            const cellShape& cell = cells[celli];

            if (!keepOrientation && !correctOrientation(points, cell))
            {
                Info<< "Inverting prism " << celli << endl;
                // Reorder prism.
                prismPoints[0] = cell[0];
                prismPoints[1] = cell[2];
                prismPoints[2] = cell[1];
                prismPoints[3] = cell[3];
                prismPoints[4] = cell[4];
                prismPoints[5] = cell[5];

                cells[celli].reset(prism, prismPoints);
            }

            celli++;

            nPrism++;
        }
        else if (elmType == MSHHEX)
        {
            storeCellInZone
            (
                regPhys,
                celli,
                physToZone,
                zoneToPhys,
                zoneCells
            );

            lineStr
                >> hexPoints[0] >> hexPoints[1]
                >> hexPoints[2] >> hexPoints[3]
                >> hexPoints[4] >> hexPoints[5]
                >> hexPoints[6] >> hexPoints[7];

            inplaceRenumber(mshToFoam, hexPoints);

            cells[celli].reset(hex, hexPoints);

            const cellShape& cell = cells[celli];

            if (!keepOrientation && !correctOrientation(points, cell))
            {
                Info<< "Inverting hex " << celli << endl;
                // Reorder hex.
                hexPoints[0] = cell[4];
                hexPoints[1] = cell[5];
                hexPoints[2] = cell[6];
                hexPoints[3] = cell[7];
                hexPoints[4] = cell[0];
                hexPoints[5] = cell[1];
                hexPoints[6] = cell[2];
                hexPoints[7] = cell[3];

                cells[celli].reset(hex, hexPoints);
            }

            celli++;

            nHex++;
        }
        else
        {
            Info<< "Unhandled element " << elmType << " at line "
                << inFile.lineNumber() << endl;
        }
    }


    inFile.getLine(line);
    IStringStream tagStr(line);
    word tag(tagStr);

    if (tag != "$ENDELM" && tag != "$EndElements")
    {
        FatalIOErrorInFunction(inFile)
            << "Did not find $ENDELM tag on line "
            << inFile.lineNumber() << exit(FatalIOError);
    }


    cells.setSize(celli);

    forAll(patchFaces, patchi)
    {
        patchFaces[patchi].shrink();
    }


    Info<< "Cells:" << endl
    << "    total:" << cells.size() << endl
    << "    hex  :" << nHex << endl
    << "    prism:" << nPrism << endl
    << "    pyr  :" << nPyr << endl
    << "    tet  :" << nTet << endl
    << endl;

    if (cells.size() == 0)
    {
        FatalIOErrorInFunction(inFile)
            << "No cells read from file " << inFile.name() << nl
            << "Does your file specify any 3D elements (hex=" << MSHHEX
            << ", prism=" << MSHPRISM << ", pyramid=" << MSHPYR
            << ", tet=" << MSHTET << ")?" << nl
            << "Perhaps you have not exported the 3D elements?"
            << exit(FatalIOError);
    }

    Info<< "CellZones:" << nl
        << "Zone\tSize" << endl;

    forAll(zoneCells, zonei)
    {
        zoneCells[zonei].shrink();

        const labelList& zCells = zoneCells[zonei];

        if (zCells.size())
        {
            Info<< "    " << zonei << '\t' << zCells.size() << endl;
        }
    }
    Info<< endl;
}

void readCells
(
    const scalar versionFormat,
    const bool keepOrientation,
    const pointField& points,
    const Map<label>& mshToFoam,
    IFstream& inFile,
    cellShapeList& cells,

    labelList& patchToPhys,
    List<DynamicList<face>>& patchFaces,

    labelList& zoneToPhys,
    List<DynamicList<label>>& zoneCells,
    Map<label> surfEntityToPhysSurface,
    Map<label> volEntityToPhysVolume
)
{
    Info<< "Starting to read cells at line " << inFile.lineNumber() << endl;

    const cellModel& hex = cellModel::ref(cellModel::HEX);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& pyr = cellModel::ref(cellModel::PYR);
    const cellModel& tet = cellModel::ref(cellModel::TET);

    face triPoints(3);
    face quadPoints(4);
    labelList tetPoints(4);
    labelList pyrPoints(5);
    labelList prismPoints(6);
    labelList hexPoints(8);


    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    label nEntities, nElems, minElemTag, maxElemTag;
    lineStr >> nEntities >> nElems >> minElemTag >> maxElemTag;

    Info<< "Cells to be read:" << nElems << endl << endl;

    // Storage for all cells. Too big. Shrink later
    cells.setSize(nElems);

    label celli = 0;
    label nTet = 0;
    label nPyr = 0;
    label nPrism = 0;
    label nHex = 0;


    // From gmsh physical region to Foam patch
    Map<label> physToPatch;

    // From gmsh physical region to Foam cellZone
    Map<label> physToZone;


    for (label entityi = 0; entityi < nEntities; entityi++)
    {
        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);

        label entityDim, entityID, regPhys, elmType, nElemInBlock, elemID;
        lineStr >> entityDim >> entityID >> elmType >> nElemInBlock;

        if (entityDim == 2)
            regPhys = surfEntityToPhysSurface[entityID];
        else if (entityDim == 3)
            regPhys = volEntityToPhysVolume[entityID];
        else
            regPhys = 0; // Points and lines don't matter to openFOAM

        // regPhys on surface elements is region number.
        if (elmType == MSHLINE)
        {
            for (label entityElm = 0; entityElm < nElemInBlock; entityElm++)
            {
                inFile.getLine(line);
                IStringStream lineStr(line);
                label meshPti;
                lineStr >> elemID >> meshPti >> meshPti;
            }
        }
        else if (elmType == MSHTRI)
        {
            for (label entityElm = 0; entityElm < nElemInBlock; entityElm++)
            {
                inFile.getLine(line);
                IStringStream lineStr(line);
                lineStr >> elemID >> triPoints[0] >> triPoints[1] >> triPoints[2];

                inplaceRenumber(mshToFoam, triPoints);

                const auto regFnd = physToPatch.cfind(regPhys);

                label patchi = -1;
                if (regFnd.good())
                {
                    // Existing patch for region
                    patchi = regFnd();
                }
                else
                {
                    // New region. Allocate patch for it.
                    patchi = patchFaces.size();

                    patchFaces.setSize(patchi + 1);
                    patchToPhys.setSize(patchi + 1);

                    Info<< "Mapping region " << regPhys << " to Foam patch "
                        << patchi << endl;
                    physToPatch.insert(regPhys, patchi);
                    patchToPhys[patchi] = regPhys;
                }

                // Add triangle to correct patchFaces.
                patchFaces[patchi].append(triPoints);
            }
        }
        else if (elmType == MSHQUAD)
        {
            for (label entityElm = 0; entityElm < nElemInBlock; entityElm++)
            {
                inFile.getLine(line);
                IStringStream lineStr(line);
                lineStr >> elemID
                    >> quadPoints[0] >> quadPoints[1] >> quadPoints[2]
                    >> quadPoints[3];

                inplaceRenumber(mshToFoam, quadPoints);

                const auto regFnd = physToPatch.cfind(regPhys);

                label patchi = -1;
                if (regFnd.good())
                {
                    // Existing patch for region
                    patchi = regFnd();
                }
                else
                {
                    // New region. Allocate patch for it.
                    patchi = patchFaces.size();

                    patchFaces.setSize(patchi + 1);
                    patchToPhys.setSize(patchi + 1);

                    Info<< "Mapping region " << regPhys << " to Foam patch "
                        << patchi << endl;
                    physToPatch.insert(regPhys, patchi);
                    patchToPhys[patchi] = regPhys;
                }

                // Add quad to correct patchFaces.
                patchFaces[patchi].append(quadPoints);
            }
        }
        else if (elmType == MSHTET)
        {
            nTet += nElemInBlock;

            for (label entityElm = 0; entityElm < nElemInBlock; entityElm++)
            {
                inFile.getLine(line);
                IStringStream lineStr(line);

                storeCellInZone
                (
                    regPhys,
                    celli,
                    physToZone,
                    zoneToPhys,
                    zoneCells
                );

                lineStr >> elemID
                    >> tetPoints[0] >> tetPoints[1] >> tetPoints[2]
                    >> tetPoints[3];

                inplaceRenumber(mshToFoam, tetPoints);

                cells[celli++].reset(tet, tetPoints);
            }
        }
        else if (elmType == MSHPYR)
        {
            nPyr += nElemInBlock;

            for (label entityElm = 0; entityElm < nElemInBlock; entityElm++)
            {
                inFile.getLine(line);
                IStringStream lineStr(line);

                storeCellInZone
                (
                    regPhys,
                    celli,
                    physToZone,
                    zoneToPhys,
                    zoneCells
                );

                lineStr >> elemID
                    >> pyrPoints[0] >> pyrPoints[1] >> pyrPoints[2]
                    >> pyrPoints[3] >> pyrPoints[4];

                inplaceRenumber(mshToFoam, pyrPoints);

                cells[celli++].reset(pyr, pyrPoints);
            }
        }
        else if (elmType == MSHPRISM)
        {
            nPrism += nElemInBlock;

            for (label entityElm = 0; entityElm < nElemInBlock; entityElm++)
            {
                inFile.getLine(line);
                IStringStream lineStr(line);

                storeCellInZone
                (
                    regPhys,
                    celli,
                    physToZone,
                    zoneToPhys,
                    zoneCells
                );

                lineStr >> elemID
                    >> prismPoints[0] >> prismPoints[1] >> prismPoints[2]
                    >> prismPoints[3] >> prismPoints[4] >> prismPoints[5];

                inplaceRenumber(mshToFoam, prismPoints);

                cells[celli].reset(prism, prismPoints);

                const cellShape& cell = cells[celli];

                if (!keepOrientation && !correctOrientation(points, cell))
                {
                    Info<< "Inverting prism " << celli << endl;
                    // Reorder prism.
                    prismPoints[0] = cell[0];
                    prismPoints[1] = cell[2];
                    prismPoints[2] = cell[1];
                    prismPoints[3] = cell[3];
                    prismPoints[4] = cell[4];
                    prismPoints[5] = cell[5];

                    cells[celli].reset(prism, prismPoints);
                }

                celli++;
            }
        }
        else if (elmType == MSHHEX)
        {
            nHex += nElemInBlock;

            for (label entityElm = 0; entityElm < nElemInBlock; entityElm++)
            {
                inFile.getLine(line);
                IStringStream lineStr(line);

                storeCellInZone
                (
                    regPhys,
                    celli,
                    physToZone,
                    zoneToPhys,
                    zoneCells
                );

                lineStr >> elemID
                    >> hexPoints[0] >> hexPoints[1]
                    >> hexPoints[2] >> hexPoints[3]
                    >> hexPoints[4] >> hexPoints[5]
                    >> hexPoints[6] >> hexPoints[7];

                inplaceRenumber(mshToFoam, hexPoints);

                cells[celli].reset(hex, hexPoints);

                const cellShape& cell = cells[celli];

                if (!keepOrientation && !correctOrientation(points, cell))
                {
                    Info<< "Inverting hex " << celli << endl;
                    // Reorder hex.
                    hexPoints[0] = cell[4];
                    hexPoints[1] = cell[5];
                    hexPoints[2] = cell[6];
                    hexPoints[3] = cell[7];
                    hexPoints[4] = cell[0];
                    hexPoints[5] = cell[1];
                    hexPoints[6] = cell[2];
                    hexPoints[7] = cell[3];

                    cells[celli].reset(hex, hexPoints);
                }

                celli++;
            }
        }
        else
        {
            Info<< "Unhandled element " << elmType << " at line "
                << inFile.lineNumber() << "in/on physical region ID: "
                << regPhys << endl;
            Info << "Perhaps you created a higher order mesh?" << endl;
        }
    }


    inFile.getLine(line);
    IStringStream tagStr(line);
    word tag(tagStr);

    if (tag != "$ENDELM" && tag != "$EndElements")
    {
        FatalIOErrorInFunction(inFile)
            << "Did not find $ENDELM tag on line "
            << inFile.lineNumber() << exit(FatalIOError);
    }


    cells.setSize(celli);

    forAll(patchFaces, patchi)
    {
        patchFaces[patchi].shrink();
    }


    Info<< "Cells:" << endl
    << "    total: " << cells.size() << endl
    << "    hex  : " << nHex << endl
    << "    prism: " << nPrism << endl
    << "    pyr  : " << nPyr << endl
    << "    tet  : " << nTet << endl
    << endl;

    if (cells.size() == 0)
    {
        FatalIOErrorInFunction(inFile)
            << "No cells read from file " << inFile.name() << nl
            << "Does your file specify any 3D elements (hex=" << MSHHEX
            << ", prism=" << MSHPRISM << ", pyramid=" << MSHPYR
            << ", tet=" << MSHTET << ")?" << nl
            << "Perhaps you have not exported the 3D elements?"
            << exit(FatalIOError);
    }

    Info<< "CellZones:" << nl
        << "Zone\tSize" << endl;

    forAll(zoneCells, zonei)
    {
        zoneCells[zonei].shrink();

        const labelList& zCells = zoneCells[zonei];

        if (zCells.size())
        {
            Info<< "    " << zonei << '\t' << zCells.size() << endl;
        }
    }
    Info<< endl;
}

void readEntities
(
    IFstream& inFile,
    Map<label>& surfEntityToPhysSurface,
    Map<label>& volEntityToPhysVolume
)
{
    label nPoints, nCurves, nSurfaces, nVolumes;
    label entityID, physicalID, nPhysicalTags;
    scalar pt; // unused scalar, gives bounding boxes of entities
    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    lineStr >> nPoints >> nCurves >> nSurfaces >> nVolumes;

    // Skip over the points, since only the full nodes list matters.
    for (label i = 0; i < nPoints; ++i)
        inFile.getLine(line);

    // Skip over the curves
    for (label i = 0; i < nCurves; ++i)
        inFile.getLine(line);

    // Read in physical surface entity groupings
    for (label i = 0; i < nSurfaces; ++i)
    {
        inFile.getLine(line);
        IStringStream lineStr(line);
        lineStr >> entityID;

        // Skip 6 useless (to us) numbers
        for (label j = 0; j < 6; ++j)
            lineStr >> pt;

        // Number of physical groups associated to this surface
        lineStr >> nPhysicalTags;
        if (nPhysicalTags > 1)
        {
            FatalIOErrorInFunction(inFile)
                << "Cannot interpret multiple physical surfaces associated"
                << " with one surface on line number " << inFile.lineNumber()
                << exit(FatalIOError);
        }

        lineStr >> physicalID;
        surfEntityToPhysSurface.insert(entityID, physicalID);
    }

    // Read in physical volume entity groupings
    for (label i = 0; i < nVolumes; ++i)
    {
        inFile.getLine(line);
        IStringStream lineStr(line);
        lineStr >> entityID;

        // Skip 6 useless (to us) numbers
        for (label j = 0; j < 6; ++j)
            lineStr >> pt;

        // Number of physical groups associated to this volume
        lineStr >> nPhysicalTags;
        if (nPhysicalTags > 1)
        {
            FatalIOErrorInFunction(inFile)
                << "Cannot interpret multiple physical volumes associated"
                << " with one volume on line number " << inFile.lineNumber()
                << exit(FatalIOError);
        }

        lineStr >> physicalID;
        volEntityToPhysVolume.insert(entityID, physicalID);
    }

    // Try to read end of section tag:
    inFile.getLine(line);
    IStringStream tagStr(line);
    word tag(tagStr);

    if (tag != "$EndEntities")
    {
        FatalIOErrorInFunction(inFile)
            << "Did not find $EndEntities tag on line "
            << inFile.lineNumber() << exit(FatalIOError);
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert a gmsh .msh file to OpenFOAM"
    );

    argList::noParallel();
    argList::addArgument(".msh file");
    argList::addBoolOption
    (
        "keepOrientation",
        "Retain raw orientation for prisms/hexs"
    );

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    // Specified region or default region
    #include "getRegionOption.H"

    if (!polyMesh::regionName(regionName).empty())
    {
        Info<< "Creating polyMesh for region " << regionName << nl;
    }

    const bool keepOrientation = args.found("keepOrientation");
    IFstream inFile(args.get<fileName>(1));

    // Storage for points
    pointField points;
    Map<label> mshToFoam;

    // Storage for all cells.
    cellShapeList cells;

    // Map from patch to gmsh physical region
    labelList patchToPhys;
    // Storage for patch faces.
    List<DynamicList<face>> patchFaces(0);

    // Map from cellZone to gmsh physical region
    labelList zoneToPhys;
    // Storage for cell zones.
    List<DynamicList<label>> zoneCells(0);

    // Name per physical region
    Map<word> physicalNames;

    // Maps from 2 and 3 dimensional entity IDs to physical region ID
    Map<label> surfEntityToPhysSurface;
    Map<label> volEntityToPhysVolume;

    // Version 1 or 2 format
    scalar versionFormat = 1;

    do
    {
        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);

        // This implies the end of while has been reached
        if (line == "")
            break;

        word tag(lineStr);

        if (tag == "$MeshFormat")
        {
            versionFormat = readMeshFormat(inFile);
        }
        else if (tag == "$PhysicalNames")
        {
            readPhysNames(inFile, physicalNames);
        }
        else if (tag == "$Entities")
        {
            // This will only happen to .msh files over version 4.
            readEntities(inFile,
                         surfEntityToPhysSurface,
                         volEntityToPhysVolume);
        }
        else if (tag == "$NOD" || tag == "$Nodes")
        {
            if (versionFormat < 4.0)
                readPointsLegacy(inFile, points, mshToFoam);
            else
                readPoints(inFile, points, mshToFoam);
        }
        else if (tag == "$ELM" || tag == "$Elements")
        {
            if (versionFormat < 4.0)
                readCellsLegacy
                (
                    versionFormat,
                    keepOrientation,
                    points,
                    mshToFoam,
                    inFile,
                    cells,
                    patchToPhys,
                    patchFaces,
                    zoneToPhys,
                    zoneCells
                );
            else
                readCells
                (
                    versionFormat,
                    keepOrientation,
                    points,
                    mshToFoam,
                    inFile,
                    cells,
                    patchToPhys,
                    patchFaces,
                    zoneToPhys,
                    zoneCells,
                    surfEntityToPhysSurface,
                    volEntityToPhysVolume
                );
        }
        else
        {
            Info<< "Skipping tag " << tag << " at line "
                << inFile.lineNumber()
                << endl;

            if (!skipSection(inFile))
            {
                break;
            }
        }
    } while (inFile.good());


    label nValidCellZones = 0;

    forAll(zoneCells, zonei)
    {
        if (zoneCells[zonei].size())
        {
            ++nValidCellZones;
        }
    }


    // Problem is that the orientation of the patchFaces does not have to
    // be consistent with the outwards orientation of the mesh faces. So
    // we have to construct the mesh in two stages:
    // 1. define mesh with all boundary faces in one patch
    // 2. use the read patchFaces to find the corresponding boundary face
    //    and repatch it.


    // Create correct number of patches
    // (but without any faces in it)
    faceListList boundaryFaces(patchFaces.size());

    wordList boundaryPatchNames(boundaryFaces.size());

    forAll(boundaryPatchNames, patchi)
    {
        boundaryPatchNames[patchi] =
            physicalNames.lookup
            (
                patchToPhys[patchi],
                polyPatch::defaultName(patchi)
            );

        Info<< "Patch " << patchi << " gets name "
            << boundaryPatchNames[patchi] << endl;
    }
    Info<< endl;

    wordList boundaryPatchTypes(boundaryFaces.size(), polyPatch::typeName);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = polyPatch::typeName;
    wordList boundaryPatchPhysicalTypes
    (
        boundaryFaces.size(),
        polyPatch::typeName
    );

    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        std::move(points),
        cells,
        boundaryFaces,
        boundaryPatchNames,
        boundaryPatchTypes,
        defaultFacesName,
        defaultFacesType,
        boundaryPatchPhysicalTypes
    );

    // Remove files now, to ensure all mesh files written are consistent.
    mesh.removeFiles();

    repatchPolyTopoChanger repatcher(mesh);

    // Now use the patchFaces to patch up the outside faces of the mesh.

    // Get the patch for all the outside faces (= default patch added as last)
    const polyPatch& pp = mesh.boundaryMesh().last();

    // Storage for faceZones.
    List<DynamicList<label>> zoneFaces(patchFaces.size());


    // Go through all the patchFaces and find corresponding face in pp.
    forAll(patchFaces, patchi)
    {
        const DynamicList<face>& pFaces = patchFaces[patchi];

        Info<< "Finding faces of patch " << patchi << endl;

        forAll(pFaces, i)
        {
            const face& f = pFaces[i];

            // Find face in pp using all vertices of f.
            label patchFacei = findFace(pp, f);

            if (patchFacei != -1)
            {
                label meshFacei = pp.start() + patchFacei;

                repatcher.changePatchID(meshFacei, patchi);
            }
            else
            {
                // Maybe internal face? If so add to faceZone with same index
                // - might be useful.
                label meshFacei = findInternalFace(mesh, f);

                if (meshFacei != -1)
                {
                    zoneFaces[patchi].append(meshFacei);
                }
                else
                {
                    WarningInFunction
                        << "Could not match gmsh face " << f
                        << " to any of the interior or exterior faces"
                        << " that share the same 0th point" << endl;
                }
            }
        }
    }
    Info<< nl;

    // Face zones
    label nValidFaceZones = 0;

    Info<< "FaceZones:" << nl
        << "Zone\tSize" << endl;

    forAll(zoneFaces, zonei)
    {
        zoneFaces[zonei].shrink();

        const labelList& zFaces = zoneFaces[zonei];

        if (zFaces.size())
        {
            ++nValidFaceZones;

            Info<< "    " << zonei << '\t' << zFaces.size() << endl;
        }
    }
    Info<< endl;


    //Get polyMesh to write to constant

    runTime.setTime(instant(0, runTime.constant()), 0);

    repatcher.repatch();

    List<cellZone*> cz;
    List<faceZone*> fz;

    // Construct and add the zones. Note that cell ordering does not change
    // because of repatch() and neither does internal faces so we can
    // use the zoneCells/zoneFaces as is.

    if (nValidCellZones > 0)
    {
        cz.setSize(nValidCellZones);

        nValidCellZones = 0;

        forAll(zoneCells, zonei)
        {
            if (zoneCells[zonei].size())
            {
                const word zoneName
                (
                    physicalNames.lookup
                    (
                        zoneToPhys[zonei],
                        "cellZone_" + Foam::name(zonei)  // default name
                    )
                );

                Info<< "Writing zone " << zonei << " to cellZone "
                    << zoneName << " and cellSet"
                    << endl;

                cellSet cset(mesh, zoneName, zoneCells[zonei]);
                cset.write();

                cz[nValidCellZones] = new cellZone
                (
                    zoneName,
                    zoneCells[zonei],
                    nValidCellZones,
                    mesh.cellZones()
                );
                nValidCellZones++;
            }
        }
    }

    if (nValidFaceZones > 0)
    {
        fz.setSize(nValidFaceZones);

        nValidFaceZones = 0;

        forAll(zoneFaces, zonei)
        {
            if (zoneFaces[zonei].size())
            {
                const word zoneName
                (
                    physicalNames.lookup
                    (
                        patchToPhys[zonei],
                        "faceZone_" + Foam::name(zonei)  // default name
                    )
                );

                Info<< "Writing zone " << zonei << " to faceZone "
                    << zoneName << " and faceSet"
                    << endl;

                faceSet fset(mesh, zoneName, zoneFaces[zonei]);
                fset.write();

                fz[nValidFaceZones] = new faceZone
                (
                    zoneName,
                    zoneFaces[zonei],
                    true, // all are flipped
                    nValidFaceZones,
                    mesh.faceZones()
                );
                nValidFaceZones++;
            }
        }
    }

    if (cz.size() || fz.size())
    {
        mesh.addZones(List<pointZone*>(0), fz, cz);
    }

    // Remove empty defaultFaces
    label defaultPatchID = mesh.boundaryMesh().findPatchID(defaultFacesName);
    if (mesh.boundaryMesh()[defaultPatchID].size() == 0)
    {
        polyPatchList newPatches((mesh.boundaryMesh().size() - 1));
        label newPatchi = 0;
        forAll(mesh.boundaryMesh(), patchi)
        {
            if (patchi != defaultPatchID)
            {
                const polyPatch& patch = mesh.boundaryMesh()[patchi];

                newPatches.set
                (
                    newPatchi,
                    patch.clone
                    (
                        mesh.boundaryMesh(),
                        newPatchi,
                        patch.size(),
                        patch.start()
                    )
                );
                ++newPatchi;
            }
        }
        repatcher.changePatches(newPatches);
    }

    // More precision (for points data)
    IOstream::minPrecision(10);

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
