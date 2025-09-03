/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "STARCDMeshWriter.H"
#include "Time.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::fileFormats::STARCDMeshWriter::findDefaultBoundary() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Find "Default_Boundary_Region" if it exists
    forAll(patches, patchi)
    {
        if (defaultBoundaryName == patches[patchi].name())
        {
            return patchi;
            break;
        }
    }

    return -1;
}


void Foam::fileFormats::STARCDMeshWriter::getCellTable()
{
    // Read constant/polyMesh/propertyName
    IOList<label> ioList
    (
        IOobject
        (
            "cellTableId",
            mesh_.time().constant(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    bool useCellZones = false;
    cellTableId_.setSize(mesh_.nCells(), -1);

    // get information from constant/polyMesh/cellTableId if possible
    if (ioList.headerOk())
    {
        if (ioList.size() == mesh_.nCells())
        {
            cellTableId_.transfer(ioList);

            if (cellTable_.empty())
            {
                Info<< "no cellTable information available" << endl;
            }
        }
        else
        {
            WarningInFunction
                << ioList.objectPath() << " has incorrect number of cells "
                << " - use cellZone information"
                << endl;

            ioList.clear();
            useCellZones = true;
        }
    }
    else
    {
        useCellZones = true;
    }


    if (useCellZones)
    {
        if (cellTable_.empty())
        {
            Info<< "created cellTable from cellZones" << endl;
            cellTable_ = mesh_;
        }

        // track if there are unzoned cells
        label nUnzoned = mesh_.nCells();

        // get the cellZone <-> cellTable correspondence
        Info<< "matching cellZones to cellTable" << endl;

        for (const cellZone& cZone : mesh_.cellZones())
        {
            if (cZone.size())
            {
                nUnzoned -= cZone.size();

                label tableId = cellTable_.findIndex(cZone.name());
                if (tableId < 0)
                {
                    dictionary dict;

                    dict.add("Label", cZone.name());
                    dict.add("MaterialType", "fluid");
                    tableId = cellTable_.append(dict);
                }

                for (const label celli : cZone)
                {
                    cellTableId_[celli] = tableId;
                }
            }
        }

        if (nUnzoned)
        {
            dictionary dict;

            dict.add("Label", "__unZonedCells__");
            dict.add("MaterialType", "fluid");
            const label tableId = cellTable_.append(dict);

            forAll(cellTableId_, i)
            {
                if (cellTableId_[i] < 0)
                {
                    cellTableId_[i] = tableId;
                }
            }
        }
    }
}


void Foam::fileFormats::STARCDMeshWriter::writeCells
(
    const fileName& prefix
) const
{
    OFstream os(starFileName(prefix, STARCDCore::CEL_FILE));
    writeHeader(os, STARCDCore::HEADER_CEL);

    //
    // Mapping between OpenFOAM and PROSTAR primitives
    //
    const Map<label> shapeLookupIndex
    {
        { cellModel::ref(cellModel::HEX).index(), STARCDCore::starcdHex },
        { cellModel::ref(cellModel::PRISM).index(), STARCDCore::starcdPrism },
        { cellModel::ref(cellModel::TET).index(), STARCDCore::starcdTet },
        { cellModel::ref(cellModel::PYR).index(), STARCDCore::starcdPyr },
    };

    const cellShapeList& shapes = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();
    const faceList& faces  = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();

    Info<< "Writing " << os.name() << " : "
        << cells.size() << " cells" << endl;

    forAll(cells, cellId)
    {
        const label tableId = cellTableId_[cellId];
        label materialType = STARCDCore::starcdFluidType; // 1(fluid)
        if (cellTable_.found(tableId))
        {
            const dictionary& dict = cellTable_[tableId];
            word matType;

            if
            (
                dict.readIfPresent("MaterialType", matType)
             && matType == "solid"
            )
            {
                materialType = STARCDCore::starcdSolidType; // 2(solid)
            }
        }

        const cellShape& shape = shapes[cellId];
        const label mapIndex = shape.model().index();

        // A registered primitive type
        if (shapeLookupIndex.found(mapIndex))
        {
            const label shapeId = shapeLookupIndex[mapIndex];
            const labelList& vrtList = shapes[cellId];

            os  << cellId + 1
                << ' ' << shapeId
                << ' ' << vrtList.size()
                << ' ' << tableId
                << ' ' << materialType;

            // Primitives have <= 8 vertices, but prevent overrun anyhow
            // indent following lines for ease of reading
            label count = 0;
            for (const label pointi : vrtList)
            {
                if ((count % 8) == 0)
                {
                    os  << nl
                        << "  " << cellId + 1;
                }
                os << ' ' << pointi + 1;
                ++count;
            }
            os << nl;
        }
        else
        {
            // Treat as general polyhedral
            const label shapeId = STARCDCore::starcdPoly;
            const labelList& cFaces  = cells[cellId];

            // create (beg,end) indices
            List<label> indices(cFaces.size() + 1);
            indices[0] = indices.size();

            label count = indices.size();
            // determine the total number of vertices
            forAll(cFaces, facei)
            {
                count += faces[cFaces[facei]].size();
                indices[facei+1] = count;
            }

            os  << cellId + 1
                << ' ' << shapeId
                << ' ' << count
                << ' ' << tableId
                << ' ' << materialType;

            // Write indices - max 8 per line
            // indent following lines for ease of reading
            count = 0;
            for (const label pointi : indices)
            {
                if ((count % 8) == 0)
                {
                    os  << nl
                        << "  " << cellId + 1;
                }
                os << ' ' << pointi;
                ++count;
            }

            // write faces - max 8 per line
            for (const label meshFace : cFaces)
            {
                face f;

                if (owner[meshFace] == cellId)
                {
                    f = faces[meshFace];
                }
                else
                {
                    f = faces[meshFace].reverseFace();
                }

                for (const label pointi : f)
                {
                    if ((count % 8) == 0)
                    {
                        os  << nl
                            << "  " << cellId + 1;
                    }

                    os << ' ' << pointi + 1;
                    ++count;
                }
            }

            os << endl;
        }
    }
}


void Foam::fileFormats::STARCDMeshWriter::writeBoundary
(
    const fileName& prefix
) const
{
    OFstream os(starFileName(prefix, STARCDCore::BND_FILE));
    writeHeader(os, STARCDCore::HEADER_BND);

    const cellShapeList& shapes = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();
    const faceList& faces  = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    //
    // Mapping between OpenFOAM and PROSTAR primitives
    // - needed for face mapping
    //
    const Map<label> shapeLookupIndex =
    {
        { cellModel::ref(cellModel::HEX).index(), STARCDCore::starcdHex },
        { cellModel::ref(cellModel::PRISM).index(), STARCDCore::starcdPrism },
        { cellModel::ref(cellModel::TET).index(), STARCDCore::starcdTet },
        { cellModel::ref(cellModel::PYR).index(), STARCDCore::starcdPyr },
    };

    Info<< "Writing " << os.name() << " : "
        << (mesh_.nFaces() - patches[0].start()) << " boundaries" << endl;


    const label defaultId = findDefaultBoundary();

    //
    // write boundary faces - skip Default_Boundary_Region entirely
    //
    label boundId = 0;
    forAll(patches, patchi)
    {
        label regionId = patchi;
        if (regionId == defaultId)
        {
            continue;  // skip - already written
        }
        else if (defaultId == -1 || regionId < defaultId)
        {
            ++regionId;
        }

        label patchStart = patches[patchi].start();
        label patchSize  = patches[patchi].size();
        word  bndType = boundaryRegion_.boundaryType(patches[patchi].name());

        for
        (
            label facei = patchStart;
            facei < (patchStart + patchSize);
            ++facei
        )
        {
            label cellId = owner[facei];
            const labelList& cFaces  = cells[cellId];
            const cellShape& shape = shapes[cellId];
            label cellFaceId = cFaces.find(facei);

            //      Info<< "cell " << cellId + 1 << " face " << facei
            //          << " == " << faces[facei]
            //          << " is index " << cellFaceId << " from " << cFaces;

            // Unfortunately, the order of faces returned by
            //   primitiveMesh::cells() is not necessarily the same
            //   as defined by primitiveMesh::cellShapes()
            // Thus, for registered primitive types, do the lookup ourselves.
            // Finally, the cellModel face number is re-mapped to the
            // STARCD local face number

            label mapIndex = shape.model().index();

            // A registered primitive type
            if (shapeLookupIndex.found(mapIndex))
            {
                const faceList sFaces = shape.faces();
                forAll(sFaces, sFacei)
                {
                    if (faces[facei] == sFaces[sFacei])
                    {
                        cellFaceId = sFacei;
                        break;
                    }
                }

                mapIndex = shapeLookupIndex[mapIndex];
                cellFaceId =
                    STARCDCore::foamToStarFaceAddr[mapIndex][cellFaceId];
            }
            // Info<< endl;

            ++boundId;

            os
                << boundId
                << ' ' << cellId + 1
                << ' ' << cellFaceId + 1
                << ' ' << regionId
                << ' ' << 0
                << ' ' << bndType.c_str()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::STARCDMeshWriter::STARCDMeshWriter
(
    const polyMesh& mesh,
    const scalar scaleFactor,
    const bool writeBndFile
)
:
    meshWriter(mesh, scaleFactor),
    writeBoundary_(writeBndFile)
{
    boundaryRegion_.readDict(mesh_);
    cellTable_.readDict(mesh_);
    getCellTable();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::STARCDMeshWriter::write(const fileName& meshName) const
{
    fileName baseName(meshName);

    if (baseName.empty())
    {
        baseName = meshWriter::defaultMeshName;

        if
        (
            mesh_.time().timeName() != "0"
         && mesh_.time().timeName() != mesh_.time().constant()
        )
        {
            baseName += "_" + mesh_.time().timeName();
        }
    }

    STARCDCore::removeFiles(baseName);

    // Points
    {
        OFstream os
        (
            starFileName(baseName, STARCDCore::VRT_FILE)
        );

        Info<< "Writing " << os.name() << " : "
            << mesh_.nPoints() << " points" << endl;

        writePoints(os, mesh_.points(), scaleFactor_);
    }

    // Cells
    writeCells(baseName);

    // Boundaries
    if (writeBoundary_)
    {
        writeBoundary(baseName);
    }

    return true;
}


// ************************************************************************* //
