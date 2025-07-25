/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    plot3dToFoam

Group
    grpMeshConversionUtilities

Description
    Plot3d mesh (ascii/formatted format) converter.

    Work in progress! Handles ascii multiblock (and optionally singleBlock)
    format.
    By default expects blanking. Use -noBlank if none.
    Use -2D \a thickness if 2D.

    Niklas Nordin has experienced a problem with lefthandedness of the blocks.
    The code should detect this automatically - see hexBlock::readPoints but
    if this goes wrong just set the blockHandedness_ variable to 'right'
    always.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IFstream.H"
#include "hexBlock.H"
#include "polyMesh.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cellShape.H"
#include "mergePoints.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Plot3d mesh (ascii/formatted format) converter"
    );
    argList::noParallel();
    argList::addArgument("PLOT3D geom file");
    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor - default is 1"
    );
    argList::addBoolOption
    (
        "noBlank",
        "Skip blank items"
    );
    argList::addBoolOption
    (
        "singleBlock",
        "Input is a single block"
    );
    argList::addOption
    (
        "2D",
        "thickness",
        "Use when converting a 2-D mesh (applied before scale)"
    );

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    const scalar scaleFactor = args.getOrDefault<scalar>("scale", 1);

    const bool readBlank = !args.found("noBlank");
    const bool singleBlock = args.found("singleBlock");
    scalar twoDThickness = -1;
    if (args.readIfPresent("2D", twoDThickness))
    {
        Info<< "Reading 2D case by extruding points by " << twoDThickness
            << " in z direction." << nl << endl;
    }


    #include "createTime.H"

    IFstream plot3dFile(args.get<fileName>(1));

    // Read the plot3d information using a fixed format reader.
    // Comments in the file are in C++ style, so the stream parser will remove
    // them with no intervention
    label nblock;

    if (singleBlock)
    {
        nblock = 1;
    }
    else
    {
        plot3dFile >> nblock;
    }

    Info<< "Reading " << nblock << " blocks" << endl;

    PtrList<hexBlock> blocks(nblock);

    {
        label nx, ny, nz;

        forAll(blocks, blockI)
        {
            if (twoDThickness > 0)
            {
                // Fake second set of points (done in readPoints below)
                plot3dFile >> nx >> ny;
                nz = 2;
            }
            else
            {
                plot3dFile >> nx >> ny >> nz;
            }

            Info<< "block " << blockI << " nx:" << nx
                << " ny:" << ny << " nz:" << nz << endl;

            blocks.set(blockI, new hexBlock(nx, ny, nz));
        }
    }

    Info<< "Reading block points" << endl;
    label sumPoints(0);
    label nMeshCells(0);

    forAll(blocks, blockI)
    {
        Info<< "block " << blockI << ":" << nl;
        blocks[blockI].readPoints(readBlank, twoDThickness, plot3dFile);
        sumPoints += blocks[blockI].nBlockPoints();
        nMeshCells += blocks[blockI].nBlockCells();
        Info<< nl;
    }

    pointField points(sumPoints);
    labelList blockOffsets(blocks.size());
    sumPoints = 0;
    forAll(blocks, blockI)
    {
        const pointField& blockPoints = blocks[blockI].points();
        blockOffsets[blockI] = sumPoints;
        forAll(blockPoints, i)
        {
            points[sumPoints++] = blockPoints[i];
        }
    }

    // From old to new master point
    labelList oldToNew;
    pointField newPoints;

    // Merge points
    mergePoints
    (
        points,
        SMALL,
        false,
        oldToNew,
        newPoints
    );

    Info<< "Merged points within " << SMALL << " distance. Merged from "
        << oldToNew.size() << " down to " << newPoints.size()
        << " points." << endl;

    // Scale the points
    if (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
    {
        newPoints *= scaleFactor;
    }

    Info<< "Creating cells" << endl;

    cellShapeList cellShapes(nMeshCells);

    const cellModel& hex = cellModel::ref(cellModel::HEX);

    label nCreatedCells = 0;

    forAll(blocks, blockI)
    {
        labelListList curBlockCells = blocks[blockI].blockCells();

        forAll(curBlockCells, blockCelli)
        {
            labelList cellPoints(curBlockCells[blockCelli].size());

            forAll(cellPoints, pointi)
            {
                cellPoints[pointi] =
                    oldToNew
                    [
                        curBlockCells[blockCelli][pointi]
                      + blockOffsets[blockI]
                    ];
            }

            // Do automatic collapse from hex.
            cellShapes[nCreatedCells].reset(hex, cellPoints, true);

            nCreatedCells++;
        }
    }

    Info<< "Creating boundary patches" << endl;

    faceListList boundary(0);
    wordList patchNames(0);
    wordList patchTypes(0);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = wallPolyPatch::typeName;
    wordList patchPhysicalTypes(0);

    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        std::move(newPoints),
        cellShapes,
        boundary,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< "Writing polyMesh" << endl;
    pShapeMesh.removeFiles();
    pShapeMesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
