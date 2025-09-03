/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    ansysToFoam

Group
    grpMeshConversionUtilities

Description
    Convert an ANSYS input mesh file (exported from I-DEAS)
    to OpenFOAM format.

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------- *\
   ------ local definitions
 \* ------------------------------------------------------------------------- */

#include <sstream>
// For EOF only
#include <cstdio>

#include "scalar.H"
#include "StringStream.H"

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "preservePatchTypes.H"
#include "cellShape.H"
#include "SLList.H"
#include "SLPtrList.H"

// Flex may use register, which is deprecated and incompatible with C++17
#pragma clang diagnostic ignored "-Wdeprecated-register"

using namespace Foam;

SLList<point> slPoints;
SLList<label> slPointMap;
label maxNodei = 0;

SLPtrList<labelList> slCellLabels;
SLList<label> slCellMap;
SLList<label> slCellType;
label maxCelli = 0;

PtrList<SLList<label>> slPatchCells;
PtrList<SLList<label>> slPatchCellFaces;

// Cell types
Map<word> cellTypes;
label currentTypei = -1;

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
#if YY_FLEX_MAJOR_VERSION <= 2 && YY_FLEX_MINOR_VERSION <= 5 && YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+
cspace                ","{space}

alpha                 [_A-Za-z]
digit                 [0-9]

identifier            {alpha}({alpha}|{digit})*
integer               {digit}+
label                 [1-9]{digit}*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))

x                     {floatNum}
y                     {floatNum}
z                     {floatNum}
value                 {floatNum}

node                  ^{space}"N"{cspace}
element               ^{space}"EN"{cspace}
bface                 ^{space}"SFE"{cspace}
elementTypeName       ^{space}"ET"{cspace}
elementType           ^{space}"TYPE"{cspace}


%%

%{
    labelList labels(8);
%}


 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

{node}{label}{cspace}{x}{cspace}{y}{cspace}{z}{space}\n {
        IStringStream nodeStream(YYText());
        char tag, c;
        label nodei;
        point node;
        nodeStream
            >> tag
            >> c >> nodei
            >> c >> node.x()
            >> c >> node.y()
            >> c >> node.z();

        if (nodei > maxNodei) maxNodei = nodei;

        slPointMap.append(nodei);
        slPoints.append(node);
    }


{element}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{space}\n {
        IStringStream elementStream(YYText());
        char tag, c;
        label celli;
        elementStream
            >> tag >> tag
            >> c >> celli
            >> c >> labels[0]
            >> c >> labels[1]
            >> c >> labels[2]
            >> c >> labels[3]
            >> c >> labels[4]
            >> c >> labels[5]
            >> c >> labels[6]
            >> c >> labels[7];

        if (celli > maxCelli) maxCelli = celli;

        slCellMap.append(celli);
        slCellLabels.append(new labelList(labels));
        slCellType.append(currentTypei);
    }


{bface}{label}{cspace}{label}{cspace}{identifier}{cspace}{integer}{cspace}{value}{space}\n {
        IStringStream bfaceStream(YYText());
        char tag, c;
        label elementi;
        label facei;
        scalar indexValue, unknown;
        bfaceStream
            >> tag >> tag >> tag
            >> c >> elementi
            >> c >> facei
            >> c >> tag >> tag >> tag >> tag
            >> c >> unknown
            >> c >> indexValue;

        label patchi = label(indexValue);

        if (patchi > slPatchCells.size())
        {
            slPatchCells.setSize(patchi);

            forAll(slPatchCells, i)
            {
                if (!slPatchCells.set(i))
                {
                    slPatchCells.set(i, new SLList<label>);
                }
            }
        }

        if (patchi > slPatchCellFaces.size())
        {
            slPatchCellFaces.setSize(patchi);

            forAll(slPatchCells, i)
            {
                if (!slPatchCellFaces.set(i))
                {
                    slPatchCellFaces.set(i, new SLList<label>);
                }
            }
        }

        slPatchCells[patchi-1].append(elementi);
        slPatchCellFaces[patchi-1].append(facei);
    }


{elementTypeName}{label}{cspace}{identifier}{space}\n {

        IStringStream elementStream(YYText());
        char tag,c;
        label cellTypei;
        word cellTypeName;
        elementStream
            >> tag >> tag           // skip 'ET'
            >> c >> cellTypei
            >> c >> cellTypeName;

        Info<< "Read typeName " << cellTypeName
            << " for type " << cellTypei << endl;

        cellTypes.insert(cellTypei, cellTypeName);
    }


{elementType}{label}{space}\n {
        IStringStream elementStream(YYText());
        char tag,c;
        label cellTypei;
        elementStream
            >> tag >> tag >> tag >> tag     // skip 'TYPE'
            >> c >> cellTypei;

        currentTypei = cellTypei;
    }



 /* ------------------------------------------------------------------------- *\
    ------ Ignore remaining space and \n s.  Any other characters are errors.
 \* ------------------------------------------------------------------------- */

.|\n {}


 /* ------------------------------------------------------------------------- *\
    ------ On EOF return to previous file, if none exists terminate.
 \* ------------------------------------------------------------------------- */

<<EOF>> {
            yyterminate();
    }
%%


#include "fileName.H"
#include <fstream>


label findFace(const polyMesh& mesh, const face& f)
{
    const labelList& pFaces = mesh.pointFaces()[f[0]];

    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        if (mesh.faces()[facei] == f)
        {
            return facei;
        }
    }

    FatalErrorInFunction
        << "Cannot find a face matching " << f
        << exit(FatalError);

    return -1;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert an ANSYS input mesh file (exported from I-DEAS)"
        " to OpenFOAM format"
    );

    argList::noParallel();
    argList::addArgument("ANSYS input file");
    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor - default is 1"
    );

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    const scalar scaleFactor = args.getOrDefault<scalar>("scale", 1);

    #include "createTime.H"

    const auto ansysFile = args.get<fileName>(1);
    std::ifstream ansysStream(ansysFile);

    if (!ansysStream)
    {
        FatalErrorInFunction
            << args.executable()
            << ": file " << ansysFile << " not found"
            << exit(FatalError);
    }

    yyFlexLexer lexer(&ansysStream);
    while (lexer.yylex() != 0)
    {}

    Info<< "Creating points" << endl;

    pointField points(slPoints.size());

    label i = 0;
    for (const point& pt : slPoints)
    {
        // Scale points for the given scale factor
        points[i++] = scaleFactor * pt;
    }


    labelList pointMap(maxNodei+1);

    i = 0;
    for (const label pointi : slPointMap)
    {
        pointMap[pointi] = i++;
    }

    Info<< "Creating cells" << endl;

    labelList cellMap(maxCelli+1);

    i = 0;
    for (const label celli : slCellMap)
    {
        cellMap[celli] = i++;
    }


    const cellModel& hex = cellModel::ref(cellModel::HEX);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& pyr = cellModel::ref(cellModel::PYR);
    const cellModel& tet = cellModel::ref(cellModel::TET);

    labelList labelsHex(8);
    labelList labelsPrism(6);
    labelList labelsPyramid(5);
    labelList labelsTet(4);

    cellShapeList cellShapes(slCellLabels.size());
    label nCells = 0;

    for (const labelList& labels : slCellLabels)
    {
        if      // Tetrahedron
        (
            labels[2] == labels[3]
         && labels[4] == labels[5]
         && labels[5] == labels[6]
         && labels[6] == labels[7]
        )
        {
            labelsTet[0] = pointMap[labels[0]];
            labelsTet[1] = pointMap[labels[1]];
            labelsTet[2] = pointMap[labels[2]];
            labelsTet[3] = pointMap[labels[4]];

            cellShapes[nCells++].reset(tet, labelsTet);
        }

        else if // Square-based pyramid
        (
            labels[4] == labels[5]
         && labels[5] == labels[6]
         && labels[6] == labels[7]
        )
        {
            labelsPyramid[0] = pointMap[labels[0]];
            labelsPyramid[1] = pointMap[labels[1]];
            labelsPyramid[2] = pointMap[labels[2]];
            labelsPyramid[3] = pointMap[labels[3]];
            labelsPyramid[4] = pointMap[labels[4]];

            cellShapes[nCells++].reset(pyr, labelsPyramid);
        }

        else if // Triangular prism
        (
            labels[2] == labels[3]
         && labels[6] == labels[7]
        )
        {
            labelsPrism[0] = pointMap[labels[0]];
            labelsPrism[1] = pointMap[labels[1]];
            labelsPrism[2] = pointMap[labels[2]];
            labelsPrism[3] = pointMap[labels[4]];
            labelsPrism[4] = pointMap[labels[5]];
            labelsPrism[5] = pointMap[labels[6]];

            cellShapes[nCells++].reset(prism, labelsPrism);
        }

        else // Hex
        {
            labelsHex[0] = pointMap[labels[0]];
            labelsHex[1] = pointMap[labels[1]];
            labelsHex[2] = pointMap[labels[2]];
            labelsHex[3] = pointMap[labels[3]];
            labelsHex[4] = pointMap[labels[4]];
            labelsHex[5] = pointMap[labels[5]];
            labelsHex[6] = pointMap[labels[6]];
            labelsHex[7] = pointMap[labels[7]];

            cellShapes[nCells++].reset(hex, labelsHex);
        }
    }


    const word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;

    // Create dummy mesh just to find out what are internal/external faces
    auto dummyMesh = autoPtr<polyMesh>::New
    (
        IOobject
        (
            "dummyMesh",
            runTime.constant(),
            runTime
        ),
        pointField(points),  // copy
        cellShapes,
        faceListList(),
        wordList(),
        wordList(),
        defaultFacesName,
        defaultFacesType,
        wordList()
    );


    // Warning: tet face order has changed between version 1.9.6 and 2.0
    //
    const label faceIndex[7][6] =
    {
        {-1, -1, -1, -1, -1, -1}, // 0
        {-1, -1, -1, -1, -1, -1}, // 1
        {-1, -1, -1, -1, -1, -1}, // 2
        {-1, -1, -1, -1, -1, -1}, // 3
        { 3,  2,  0, -1,  1, -1}, // tet (version 2.0)
        { 0,  4,  3, -1,  2,  1}, // prism
        { 4,  2,  1,  3,  0,  5}, // hex
    };

    Info<< "Creating boundary patches" << endl;

    faceListList boundary(slPatchCells.size());
    wordList patchNames(slPatchCells.size());

    forAll(slPatchCells, patchi)
    {
        SLList<face> patchFaces;

        auto cellIter = slPatchCells[patchi].cbegin();
        auto faceIter = slPatchCellFaces[patchi].cbegin();

        for
        (
            ;
            cellIter.good() && faceIter.good();
            ++cellIter, ++faceIter
        )
        {
            const cellShape& shape = cellShapes[cellMap[cellIter()]];

            patchFaces.append
            (
                shape.faces()
                [
                    faceIndex
                        [shape.nFaces()]
                        [faceIter()-1]
                ]
            );
        }

        boundary[patchi] = patchFaces;
        patchNames[patchi] = polyPatch::defaultName(patchi + 1);
    }


    //
    // Lookup the face labels for all the boundary faces
    //
    labelListList boundaryFaceLabels(boundary.size());
    forAll(boundary, patchi)
    {
        const faceList& bFaces = boundary[patchi];
        labelList& bFaceLabels = boundaryFaceLabels[patchi];
        bFaceLabels.setSize(bFaces.size());
        forAll(bFaces, i)
        {
            bFaceLabels[i] = findFace(dummyMesh(), bFaces[i]);
        }
    }


    // Now split the boundary faces into external and internal faces. All
    // faces go into faceZones and external faces go into patches.
    List<faceList> patchFaces(slPatchCells.size());
    labelList patchNFaces(slPatchCells.size(), Zero);
    forAll(boundary, patchi)
    {
        const faceList& bFaces = boundary[patchi];
        const labelList& bFaceLabels = boundaryFaceLabels[patchi];

        patchFaces[patchi].setSize(bFaces.size());

        forAll(bFaces, i)
        {
            if (!dummyMesh().isInternalFace(bFaceLabels[i]))
            {
                patchFaces[patchi][patchNFaces[patchi]++] = bFaces[i];
            }
        }
        patchFaces[patchi].setSize(patchNFaces[patchi]);

        Info<< "Patch " << patchi << " named " << patchNames[patchi]
            << ": " << boundary[patchi].size() << " faces" << endl;
    }

    // We no longer need the dummyMesh
    dummyMesh.clear();


    Info<< "ansysToFoam: " << endl
        << "Ansys file format does not provide information about the type of "
        << "the patch (eg. wall, symmetry plane, cyclic etc)." << endl
        << "All the patches have been created "
        << "as type patch. Please reset after mesh conversion as necessary."
        << endl;

    PtrList<dictionary> patchDicts;

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );

    // Add information to dictionary
    forAll(patchNames, patchi)
    {
        if (!patchDicts.set(patchi))
        {
            patchDicts.set(patchi, new dictionary());
        }
        // Add but not overwrite
        patchDicts[patchi].add("type", polyPatch::typeName, false);
    }


    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(points),
        cellShapes,
        patchFaces,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );


    if (cellTypes.size() || patchNames.size())
    {
        DynamicList<pointZone*> pz;
        DynamicList<faceZone*> fz;
        DynamicList<cellZone*> cz;

        // FaceZones
        forAll(boundaryFaceLabels, patchi)
        {
            if (boundaryFaceLabels[patchi].size())
            {
                // Re-do the boundaryFaceLabels since the boundary face
                // labels will be different on the pShapeMesh.
                const faceList& bFaces = boundary[patchi];
                labelList& bFaceLabels = boundaryFaceLabels[patchi];
                forAll(bFaceLabels, i)
                {
                    bFaceLabels[i] = findFace(pShapeMesh, bFaces[i]);
                }

                Info<< "Creating faceZone " <<  patchNames[patchi]
                    << " with " << bFaceLabels.size() << " faces" << endl;

                fz.append
                (
                    new faceZone
                    (
                        patchNames[patchi],
                        bFaceLabels,
                        false, // none are flipped
                        fz.size(),
                        pShapeMesh.faceZones()
                    )
                );
            }
        }


        // CellZones
        labelList types = cellTypes.sortedToc();

        forAll(types, typei)
        {
            const label cellType = types[typei];

            // Pick up cells in zone
            DynamicList<label> addr;

            auto cellMapIter = slCellMap.cbegin();
            auto typeIter = slCellType.cbegin();

            for
            (
                ;
                typeIter.good();
                ++typeIter, ++cellMapIter
            )
            {
                if (typeIter() == cellType)
                {
                    addr.append(cellMap[cellMapIter()]);
                }
            }

            Info<< "Creating cellZone " << cellTypes[cellType]
                << " with " << addr.size() << " cells" << endl;

            cz.append
            (
                new cellZone
                (
                    cellTypes[cellType],
                    addr,
                    typei,
                    pShapeMesh.cellZones()
                )
            );
        }

        pShapeMesh.addZones(pz, fz, cz);
    }


    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< "Writing polyMesh" << endl;
    pShapeMesh.removeFiles();
    pShapeMesh.write();

    Info<< nl << "End" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of ansysToFoam.L
 \* ------------------------------------------------------------------------- */
