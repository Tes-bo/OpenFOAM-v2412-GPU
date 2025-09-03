/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    blockMesh

Group
    grpMeshGenerationUtilities

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
      - \c system/blockMeshDict
      - \c system/\<region\>/blockMeshDict
      - \c constant/polyMesh/blockMeshDict
      - \c constant/\<region\>/polyMesh/blockMeshDict

Usage
    \b blockMesh [OPTION]

    Options:
      - \par -write-obj
        Write topology as a set of edges in OBJ format and exit.

      - \par -write-vtk
        Write topology as VTK file (xml, ascii) and exit.

      - \par -merge-points
        Merge points instead of default topological merge

      - \par -region \<name\>
        Specify alternative mesh region.

      - \par -dict \<filename\>
        Alternative dictionary for the block mesh description.

      - \par -sets
        Write cellZones as cellSets too (for processing purposes)

      - \par -no-clean
        Do not remove polyMesh/ directory or files

      - \par -time
        Write resulting mesh to a time directory (instead of constant)

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "foamVtkInternalMeshWriter.H"
#include "foamVtkSurfaceWriter.H"
#include "attachPolyTopoChanger.H"
#include "polyTopoChange.H"
#include "cyclicPolyPatch.H"
#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "wordPair.H"
#include "slidingInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Block mesh generator.\n"
        " \n"
        "  The ordering of vertex and face labels within a block as shown "
        "below.\n"
        "  For the local vertex numbering in the sequence 0 to 7:\n"
        "    Faces 0, 1 (x-direction) are left, right.\n"
        "    Faces 2, 3 (y-direction) are front, back.\n"
        "    Faces 4, 5 (z-direction) are bottom, top.\n"
        " \n"
        "                        7 ---- 6\n"
        "                 f5     |\\     :\\     f3\n"
        "                 |      | 4 ---- 5     \\\n"
        "                 |      3.|....2 |      \\\n"
        "                 |       \\|     \\|      f2\n"
        "                 f4       0 ---- 1\n"
        "    Y  Z\n"
        "     \\ |                f0 ------ f1\n"
        "      \\|\n"
        "       o--- X\n"
    );

    argList::noParallel();
    argList::noFunctionObjects();

    argList::addBoolOption
    (
        "write-obj",
        "Write block edges and centres as obj files and exit",
        true  // (old) mark as advanced option. -write-vtk is preferred
    );
    argList::addOptionCompat("write-obj", {"blockTopology", 1912});

    argList::addBoolOption
    (
        "write-vtk",
        "Write topology as VTU file and exit"
    );
    argList::addVerboseOption
    (
        "Force verbose output. (Can be used multiple times)"
    );

    argList::addBoolOption
    (
        "merge-points",
        "Geometric point merging instead of topological merging"
        " [default for 1912 and earlier]."
        // NOTE: " Slower, fails with high-aspect cells."
    );
    argList::addBoolOption
    (
        "no-clean",
        "Do not remove polyMesh/ directory or files"
    );
    argList::addOptionCompat("no-clean", {"noClean", -2006});

    argList::addOption("dict", "file", "Alternative blockMeshDict");
    argList::addBoolOption
    (
        "sets",
        "Write cellZones as cellSets too (for processing purposes)"
    );
    argList::addOption
    (
        "time",
        "time",
        "Specify a time to write mesh to (default: constant)"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Remove old files, unless disabled
    const bool removeOldFiles = !args.found("no-clean");

    // Write cellSets
    const bool writeCellSets = args.found("sets");

    // Default merge (topology), unless otherwise specified
    const blockMesh::mergeStrategy strategy =
    (
        args.found("merge-points")
      ? blockMesh::MERGE_POINTS
      : blockMesh::DEFAULT_MERGE
    );

    // Specified region or default region
    #include "getRegionOption.H"

    if (!polyMesh::regionName(regionName).empty())
    {
        Info<< nl << "Generating mesh for region " << regionName << nl;
    }

    // Instance for resulting mesh
    bool useTime = false;
    word meshInstance(runTime.constant());

    if
    (
        args.readIfPresent("time", meshInstance)
     && runTime.constant() != meshInstance
    )
    {
        // Verify that the value is actually good
        scalar timeValue;

        useTime = readScalar(meshInstance, timeValue);
        if (!useTime)
        {
            FatalErrorInFunction
                << "Bad input value: " << meshInstance
                << "Should be a scalar or 'constant'"
                << nl << endl
                << exit(FatalError);
        }
    }


    // Locate appropriate blockMeshDict
    #include "findBlockMeshDict.H"

    blockMesh blocks(meshDict, regionName, strategy, args.verbose());

    if (!blocks.good())
    {
        // Could/should be Fatal?

        WarningIn(args.executable())
            << "Did not generate any blocks. Stopping." << nl << endl;

        return 1;
    }


    bool quickExit = false;

    if (args.found("write-obj"))
    {
        quickExit = true;
        Info<< nl;
        #include "blockMeshOBJ.H"
    }

    if (args.found("write-vtk"))
    {
        quickExit = true;
        Info<< nl;
        #include "blockMeshVTK.H"
    }

    if (quickExit)
    {
        Info<< "\nEnd\n" << endl;
        return 0;
    }


    // Instance for resulting mesh
    if (useTime)
    {
        Info<< "Writing polyMesh to " << meshInstance << nl << endl;

        // Make sure that the time is seen to be the current time.
        // This is the logic inside regIOobject that resets the instance
        // to the current time before writing
        runTime.setTime(instant(meshInstance), 0);
    }

    if (removeOldFiles)
    {
        #include "cleanMeshDirectory.H"
    }


    // Ensure we get information messages, even if turned off in dictionary
    blocks.verbose(true);

    autoPtr<polyMesh> meshPtr =
        blocks.mesh
        (
            IOobject(regionName, meshInstance, runTime)
        );

    polyMesh& mesh = *meshPtr;

    // Merge patch pairs (dictionary entry "mergePatchPairs")
    #include "mergePatchPairs.H"

    // Handle cyclic patches
    #include "handleCyclicPatches.H"

    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< nl << "Writing polyMesh with "
        << mesh.cellZones().size() << " cellZones";

    if (writeCellSets && !mesh.cellZones().empty())
    {
        Info<< " (written as cellSets too)";
    }
    Info<< endl;

    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    if (writeCellSets)
    {
        for (const cellZone& cz : mesh.cellZones())
        {
            cellSet(mesh, cz.name(), cz).write();
        }
    }

    #include "printMeshSummary.H"

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
