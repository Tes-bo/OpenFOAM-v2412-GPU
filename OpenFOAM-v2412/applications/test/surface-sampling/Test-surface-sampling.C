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
    Test-surface-sampling

Description
    Simple test of surface sampling, including timings

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "clockTime.H"

#include "fileName.H"
#include "sampledSurfaces.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "profiling.H"
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
    #include "addProfilingOption.H"

    argList::addOption
    (
        "sample",
        "file",
        "Name of surface sampling to use (default: test-sample)"
    );
    argList::addOption("output", "Begin output iteration (default: 10000)");
    argList::addOption("count", "Number of writes (default: 1)");

    #include "setRootCase.H"
    #include "createTime.H"

    const int verbose = args.verbose();
    const label firstOutput = args.getOrDefault("output", 10000);
    const label nOutput = args.getOrDefault("count", 1);

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
            IOobject::MUST_READ
        ),
        false
    );
    mesh.init(true);   // initialise all (lower levels and current)
    Info<< endl;

    // Like "setSystemMeshDictionaryIO.H"
    IOobject dictIO = IOobject::selectIO
    (
        IOobject
        (
            "test-sample",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        args.getOrDefault<fileName>("sample", "")
    );

    dictionary dictContents = IOdictionary(dictIO);
    const dictionary* sampleDict = nullptr;

    if (!dictContents.empty())
    {
        // Either have 'regular form' (from functionObjects)
        // in which the sample surfaces are buried one layer deep
        // or a flattened dictionary

        if (dictContents.front()->dictPtr())
        {
            // Appears to be a dictionary of contents
            // - get the first sub-dictionary with the correct "type"

            for (const entry& e : dictContents)
            {
                const dictionary* dptr = e.dictPtr();

                if
                (
                    dptr
                 &&
                    (
                        sampledSurfaces::typeName
                     == dptr->getOrDefault<word>("type", word::null)
                    )
                )
                {
                    sampleDict = dptr;
                    break;
                }
            }
        }
        else
        {
            // Probably a flattened dictionary,
            // just check directly

            if
            (
                sampledSurfaces::typeName
             == dictContents.getOrDefault<word>("type", word::null)
            )
            {
                sampleDict = &dictContents;
            }
        }
    }

    if (!sampleDict)
    {
        FatalErrorInFunction
            << "Dictionary does not appear to contain type:"
            << sampledSurfaces::typeName << nl
            << "    " << dictIO.objectRelPath() << nl
            << exit(FatalError);
    }

    // Construct from Time and dictionary, without loadFromFiles
    sampledSurfaces sampling("test-sample", runTime, *sampleDict);

    #if (OPENFOAM <= 2306)
    Info<< "Loaded " << sampling.size() << " surface samplers, fields: "
        << flatOutput(sampleDict->getOrDefault<wordRes>("fields", wordRes()))
        << nl;
    #else
    Info<< "Loaded " << sampling.size() << " surface samplers, fields: "
        << flatOutput(sampling.fieldNames())
        << nl;
    #endif

    if (sampling.empty())
    {
        FatalErrorInFunction
            << "No surface samplers loaded" << nl
            << "    " << dictIO.objectRelPath() << nl
            << exit(FatalError);
    }


    // Manually read and load files

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read GeometricFields
    Info<< nl << "Load fields" << nl;

    #if (OPENFOAM <= 2306)
    // List of stored objects to clear after (as required)
    wordHashSet allFields(objects.names());
    LIFOStack<regIOobject*> storedObjects;

    #undef  ReadFields
    #define ReadFields(FieldType)                                             \
    readFields<FieldType>(mesh, objects, allFields, storedObjects);
    #else
    // List of stored objects to clear after (as required)
    DynamicList<regIOobject*> storedObjects;

    #undef  ReadFields
    #define ReadFields(FieldType)                                             \
    readFields<FieldType>(mesh, objects, predicates::always{}, storedObjects);
    #endif

    // Read volFields
    ReadFields(volScalarField);
    ReadFields(volVectorField);
    ReadFields(volSphericalTensorField);
    ReadFields(volSymmTensorField);
    ReadFields(volTensorField);

    // Set fields to AUTO_WRITE (not really necessary for sampling...)
    for (regIOobject* io : storedObjects)
    {
        io->writeOpt(IOobjectOption::AUTO_WRITE);
    }

    Info<< nl
        << "Start " << nOutput << " times starting at "
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
        sampling.write();
    }

    if (verbose) Info<< nl;
    Info<< nl << "Writing took "
        << timing.timeIncrement() << "s" << endl;

    //TBD    profiling::writeNow();

    // Info<< nl
    //     << "Cleanup newly generated files with" << nl << nl
    //     << "    foamListTimes -rm -time "
    //     << firstOutput << ":" << nl
    //     << "    foamListTimes -processor -rm -time "
    //     << firstOutput << ":" << nl;


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
