/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-fileHandler-writing

Description
    Simple test of file writing, including timings

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "clockTime.H"

#include "fileName.H"
#include "fileOperation.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "ReadFields.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Loads mesh and fields from latest time and writes multiple times"
    );

    argList::noFunctionObjects();  // Disallow function objects
    argList::addVerboseOption("additional verbosity");
    argList::addOption
    (
        "output",
        "N",
        "Begin output iteration (default: 10000)"
    );
    argList::addOption
    (
        "count",
        "N",
        "Number of writes (default: 1)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const label firstOutput = args.getOrDefault("output", 10000);
    const label nOutput = args.getOrDefault("count", 1);

    const int verbose = args.verbose();

    // Select latestTime, including 0 and constant
    {
        const auto& times = runTime.times();
        const label timeIndex = (times.size()-1);

        if (timeIndex < 0)
        {
            FatalErrorInFunction
                << "No times!"
                << exit(FatalError);
        }

        runTime.setTime(times[timeIndex], timeIndex);
    }

    // #include "createMesh.H"

    Info << "Create mesh time = " << runTime.timeName() << nl;

    fvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        false
    );
    mesh.init(true);   // initialise all (lower levels and current)
    Info<< endl;


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // List of stored objects to clear after (as required)
    DynamicList<regIOobject*> storedObjects;


    // Read GeometricFields
    Info<< nl << "Load fields" << nl;

    #undef  ReadFields
    #define ReadFields(FieldType)                                             \
    readFields<FieldType>(mesh, objects, predicates::always{}, storedObjects);

    // Read volFields
    ReadFields(volScalarField);
    ReadFields(volVectorField);
    ReadFields(volSphericalTensorField);
    ReadFields(volSymmTensorField);
    ReadFields(volTensorField);

    // Set fields to AUTO_WRITE
    for (regIOobject* io : storedObjects)
    {
        io->writeOpt(IOobjectOption::AUTO_WRITE);
    }

    Info<< nl
        << "Writing " << nOutput << " times starting at "
        << firstOutput << nl;

    clockTime timing;

    if (verbose) Info<< "Time:";

    for
    (
        label timeIndex = firstOutput, count = 0;
        count < nOutput;
        ++timeIndex, ++count
    )
    {
        runTime.setTime(timeIndex, timeIndex);
        if (verbose) Info<< ' ' << runTime.timeName() << flush;
        runTime.writeNow();
    }

    if (verbose) Info<< nl;
    Info<< nl << "Writing took "
        << timing.timeIncrement() << "s" << endl;


    Info<< nl
        << "Cleanup newly generated files with" << nl << nl
        << "    foamListTimes -rm -time "
        << firstOutput << ":" << nl
        << "    foamListTimes -processor -rm -time "
        << firstOutput << ":" << nl;


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
