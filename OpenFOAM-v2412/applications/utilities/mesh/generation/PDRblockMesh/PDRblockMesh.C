/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    PDRblockMesh

Group
    grpMeshGenerationUtilities

Description
    A specialized single-block mesh generator for a rectilinear mesh
    in x-y-z.

    Uses the mesh description found in
      - \c system/PDRblockMeshDict

Usage
    \b PDRblockMesh [OPTION]

    Options:
      - \par -dict \<filename\>
        Alternative dictionary for the mesh description.

      - \par -no-clean
        Do not remove polyMesh/ directory or files

      - \par -time
        Write resulting mesh to a time directory (instead of constant)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "PDRblock.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OSspecific.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A block mesh generator for a rectilinear mesh in x-y-z.\n"
        "  The ordering of vertex and face labels within a block as shown "
        "below.\n"
        "  For the local vertex numbering in the sequence 0 to 7:\n"
        "    Faces 0, 1  ==  x-min, x-max.\n"
        "    Faces 2, 3  ==  y-min, y-max.\n"
        "    Faces 4, 5  ==  z-min, z-max.\n"
        "\n"
        "                        7 ---- 6\n"
        "                 f5     |\\     |\\     f3\n"
        "                 |      | 4 ---- 5     \\\n"
        "                 |      3 |--- 2 |      \\\n"
        "                 |       \\|     \\|      f2\n"
        "                 f4       0 ---- 1\n"
        "    Y  Z\n"
        "     \\ |                f0 ------ f1\n"
        "      \\|\n"
        "       O--- X\n"
    );

    argList::noParallel();
    argList::noFunctionObjects();

    argList::addBoolOption
    (
        "no-clean",
        "Do not remove polyMesh/ directory or files"
    );
    argList::addOptionCompat("no-clean", {"noClean", -2006});

    argList::addBoolOption
    (
        "no-outer",
        "Create without any other region"
    );
    argList::addBoolOption
    (
        "print-dict",
        "Print blockMeshDict equivalent and exit"
    );
    argList::addBoolOption
    (
        "write-dict",
        "Write system/blockMeshDict.PDRblockMesh and exit"
    );

    argList::addOption("dict", "file", "Alternative PDRblockMeshDict");
    argList::addOption
    (
        "time",
        "time",
        "Specify a time to write mesh to (default: constant)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    // Remove old files, unless disabled
    const bool removeOldFiles = !args.found("no-clean");

    // Suppress creation of the outer region
    const bool noOuterRegion = args.found("no-outer");

    const word regionName(polyMesh::defaultRegion);


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


    // Locate appropriate PDRblockMeshDict
    const word dictName("PDRblockMeshDict");
    #include "setSystemRunTimeDictionaryIO.H"

    IOdictionary meshDict(dictIO);

    Info<< "Creating PDRblockMesh from "
        << dictIO.objectRelPath() << endl;

    // Always start from a PDRblock
    PDRblock blkMesh(meshDict, true);

    if (args.found("print-dict"))
    {
        Info<< nl << "Equivalent blockMeshDict" << nl << nl;

        blkMesh.blockMeshDict(Info, true);

        Info<< "\nEnd\n" << endl;
        return 0;
    }

    if (args.found("write-dict"))
    {
        // Generate system/blockMeshDict and exit
        blkMesh.writeBlockMeshDict
        (
            IOobject
            (
                "blockMeshDict.PDRblockMesh",
                runTime.system(),   // instance
                runTime,            // registry
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

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


    Info<< nl << "Creating polyMesh from PDRblockMesh" << endl;
    if (noOuterRegion)
    {
        Info<< "Outer region disabled, using ijk generation" << nl;
    }

    autoPtr<polyMesh> meshPtr =
    (
        args.found("no-outer")
      ? blkMesh.innerMesh(IOobject(regionName, meshInstance, runTime))
      : blkMesh.mesh(IOobject(regionName, meshInstance, runTime))
    );

    polyMesh& mesh = *meshPtr;

    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< nl << "Writing polyMesh with "
        << mesh.cellZones().size() << " cellZones" << endl;

    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    #include "printMeshSummary.H"

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
