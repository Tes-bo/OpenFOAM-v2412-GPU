/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Applications
    PDRsetFields

Description
    Preparation of fields for PDRFoam

SourceFiles
    PDRsetFields.C

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"

#include "PDRsetFields.H"
#include "PDRlegacy.H"
#include "PDRutils.H"
#include "IOmanip.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Processes a set of geometrical obstructions to determine the"
        " equivalent blockage effects when setting cases for PDRFoam"
    );
    argList::noParallel();
    argList::noFunctionObjects();

    argList::addOption
    (
        "time",
        "time",
        "Specify a time"
    );

    argList::addOption("dict", "file", "Alternative PDRsetFieldsDict");

    argList::addBoolOption
    (
        "legacy",
        "Force use of legacy obstacles table"
    );

    argList::addDryRunOption
    (
        "Read obstacles and write VTK only"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("PDRsetFieldsDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary setFieldsDict(dictIO);

    const fileName& casepath = runTime.globalPath();

    pars.timeName = "0";
    args.readIfPresent("time", pars.timeName);

    // Program parameters (globals)
    pars.read(setFieldsDict);

    if (args.found("legacy"))
    {
        pars.legacyObsSpec = true;
    }

    // Always have the following:
    // 0 = blockedFaces patch (no wall functions)
    // 1 = mergingFaces patch
    // 2 = wallFaces patch

    DynamicList<PDRpatchDef> patches;
    patches.resize(PDRpatchDef::NUM_PREDEFINED);

    for
    (
        PDRpatchDef::predefined predef :
        {
            PDRpatchDef::BLOCKED_FACE,
            PDRpatchDef::MERGING_PATCH,
            PDRpatchDef::WALL_PATCH,
        }
    )
    {
        patches[predef] = PDRpatchDef::names[predef];
    }


    // Dimensions and grid points for i-j-k domain
    PDRblock pdrBlock;

    if (pars.legacyMeshSpec)
    {
        PDRlegacy::read_mesh_spec(casepath, pdrBlock);
    }
    else
    {
        IOdictionary iodict
        (
            IOobject
            (
                "PDRblockMeshDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        pdrBlock.read(iodict);

        #ifdef FULLDEBUG
        PDRlegacy::print_info(pdrBlock);
        #endif
    }

    // Storage for obstacles and cylinder-like obstacles
    DynamicList<PDRobstacle> obstacles, cylinders;

    // Read in obstacles
    const scalar volObstacles =
    (
        pars.legacyObsSpec
      ? PDRobstacle::legacyReadFiles
        (
            pars.obsfile_dir, pars.obsfile_names,
            pdrBlock.bounds(),
            obstacles,
            cylinders
        )
      : PDRobstacle::readFiles
        (
            pars.obsfile_dir, pars.obsfile_names,
            pdrBlock.bounds(),
            obstacles,
            cylinders
        )
    );


    PDRobstacle::generateVtk(casepath/"VTK", obstacles, cylinders);

    if (args.dryRun())
    {
        Info<< nl
            << "dry-run: stopping after reading/writing obstacles" << nl
            << "\nEnd\n" << nl;
        return 0;
    }


    // Bookkeeping of the ranges within the obstacle list

    // Original blockage at the start
    const labelRange origBlocks(0, obstacles.size());

    // Intersection blockage
    labelRange interBlocks(origBlocks.end_value(), 0);

    scalar volSubtract = 0;

    // Do binary intersections between blocks and cylinders (or diag-beam)
    // by creating -ve blocks at the overlap

    labelRange int1Blocks(origBlocks.end_value(), 0);

    if (pars.overlaps % 2 > 0)
    {
        Info<< "    block/cylinder intersections" << endl;

        label nblocked = obstacles.size();

        volSubtract += block_cylinder_overlap(obstacles, origBlocks, cylinders);

        nblocked = (obstacles.size() - nblocked);

        interBlocks += nblocked;
        int1Blocks += nblocked;
    }

    // Do binary intersections between blocks
    // by creating -ve blocks at the overlap

    labelRange int2Blocks(int1Blocks.end_value(), 0);
    if (pars.overlaps % 4 > 1)
    {
        Info<< "    block/block intersections" << endl;

        label nblocked = obstacles.size();

        volSubtract += block_overlap(obstacles, origBlocks, 1.0);

        nblocked = (obstacles.size() - nblocked);

        interBlocks += nblocked;
        int2Blocks += nblocked;
    }

    // Correct for triple intersections
    // by looking for overlaps between the -ve blocks just created

    labelRange int3Blocks(int2Blocks.end_value(), 0);
    if (pars.overlaps % 8 > 3)
    {
        Info<< "    triple intersections" << endl;

        label nblocked = obstacles.size();

        volSubtract += block_overlap(obstacles, interBlocks, 1.0/3.0);

        nblocked = (obstacles.size() - nblocked);

        interBlocks += nblocked;
        int3Blocks += nblocked;
    }


    // The field arrays, in one structure pass around easily
    PDRarrays arr(pdrBlock);

    Info<< "Apply blockage" << endl;

    // Create blockage and count arrays by working through
    // real and extra blocks and cylinders

    // User-defined negative blocks. Use "sign" to distinguish
    if (origBlocks.size())
    {
        Info<< "    negative blocks: " << origBlocks.size() << nl;

        for (const PDRobstacle& obs : obstacles.slice(origBlocks))
        {
            arr.addBlockage(obs, patches, -1);
        }
    }

    // Do the intersection blocks positive and negative
    // These are done first so that negative area blockage cancels positive

    if (interBlocks.size())
    {
        Info<< "    blocks " << interBlocks.size() << nl;

        for (const PDRobstacle& obs : obstacles.slice(interBlocks))
        {
            arr.addBlockage(obs, patches, 0);
        }
    }

    // The positive real bocks
    if (origBlocks.size())
    {
        Info<< "    positive blocks: " << origBlocks.size() << nl;

        for (const PDRobstacle& obs : obstacles.slice(origBlocks))
        {
            arr.addBlockage(obs, patches, 1);
        }
    }

    // The cylinders
    if (cylinders.size())
    {
        Info<< "    cylinders: " << cylinders.size() << nl;

        for (const PDRobstacle& obs : cylinders)
        {
            arr.addCylinder(obs);
        }
    }

    // Calculation of the fields of drag, turbulence
    // generation and combustion enhancement

    arr.blockageSummary();

    // Mapping of OpenFOAM cells/faces to i-j-k indices
    PDRmeshArrays meshIdx;
    meshIdx.gridPointRelTol = pars.gridPointTol;

    meshIdx.read(runTime, pdrBlock);

    PDRarrays::calculateAndWrite(arr, meshIdx, casepath, patches);

    Info<< nl
        << setw(6) << origBlocks.size() << " blocks and "
        << cylinders.size() << " cylinders/diagonal blocks" << nl;

    Info<< setw(6) << int2Blocks.size()
        << " intersections amongst blocks" << nl;

    Info<< setw(6) << int1Blocks.size()
        << " intersections between blocks and cyl/beams" << nl;

    Info<< setw(6) << int1Blocks.size()
        << "/3 triple intersections" << nl;

    Info<< "Volume of obstacles read in: " << volObstacles
        << ", volume of intersections: " << volSubtract << nl;

    Info<< nl << "After corrections:" << nl;
    arr.blockageSummary();

    Info<< nl << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
