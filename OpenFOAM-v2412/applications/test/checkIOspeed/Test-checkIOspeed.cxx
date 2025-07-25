/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Test-checkIOspeed

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
#include "globalIndex.H"
#include "volFields.H"
#include "IOField.H"
#include "PDRblock.H"

// Not really great since CoherentMesh only works with reading!
#ifdef USE_COHERENT
#include "OFCstream.H"
#include "SliceStreamRepo.H"
#endif

#include <numeric>

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote("Rewrites fields multiple times");

    argList::noFunctionObjects();  // Disallow function objects
    argList::noCheckProcessorDirectories();

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
    argList::addOption
    (
        "fields",
        "N",
        "Number of fields to write (default: 1)"
    );
    argList::addOption
    (
        "global",
        "N",
        "Global field size"
    );
    argList::addOption
    (
        "local",
        "N",
        "Local fields size (default: 1000)"
    );
    argList::addOption
    (
        "mesh",
        "(nx ny nz)",
        "Create with a mesh"
    );
    argList::addOption
    (
        "exclude",
        "(int ... )",
        "zero-sized on ranks with specified modulo"
    );

    argList::addBoolOption("coherent", "Force coherent output");

    #include "setRootCase.H"
    #include "createTime.H"

    const label firstOutput = args.getOrDefault("output", 10000);
    const label nOutput = args.getOrDefault("count", 1);
    const label nFields = args.getOrDefault("fields", 1);
    labelVector meshCells(0, 0, 0);

    const int verbose = args.verbose();

    const bool useCoherent = args.found("coherent");

    labelList excludes;
    args.readListIfPresent("exclude", excludes);

    bool writeOnProc = true;

    const label myProci = UPstream::myProcNo();

    for (const label excl : excludes)
    {
        if (excl > 1 && myProci > 0 && (myProci % excl) == 0)
        {
            writeOnProc = false;
            break;
        }
    }

    const label nProcsEff =
        returnReduce((writeOnProc ? 1 : 0), sumOp<label>());

    Info<< "Output " << nProcsEff
        << "/" << UPstream::nProcs()
        << " ranks" << nl;


    if (args.readIfPresent("mesh", meshCells))
    {
        if (!writeOnProc)
        {
            meshCells = Zero;
        }

        PDRblock block(boundBox(point::zero, point::one), meshCells);

        // Advance time
        // - coherent currently still needs to read the mesh itself!
        runTime.setTime(firstOutput, firstOutput);

        IOobject meshIO
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        autoPtr<polyMesh> pmeshPtr(block.innerMesh(meshIO));

        fvMesh mesh
        (
            meshIO,
            pointField(pmeshPtr->points()),
            faceList(pmeshPtr->faces()),
            labelList(pmeshPtr->faceOwner()),
            labelList(pmeshPtr->faceNeighbour())
        );

        pmeshPtr.reset(nullptr);

        const label fieldSize = mesh.nCells();

        const globalIndex giCells(fieldSize);

        // Create fields
        Info<< nl << "Create " << nFields << " fields" << nl
            << "field-size:" << fieldSize
            << " total-size:" << giCells.totalSize() << nl;

        // Dimensioned field (no proc boundaries)

        PtrList<volScalarField::Internal> fields(nFields);

        {
            IOobject io
            (
                "field",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            forAll(fields, fieldi)
            {
                io.resetHeader("field" + Foam::name(fieldi));

                fields.set
                (
                    fieldi,
                    new volScalarField::Internal(io, mesh, dimless)
                );

                auto& fld = fields[fieldi];
                std::iota(fld.begin(), fld.end(), scalar(0));
            }
        }

        IOstreamOption streamOpt(IOstreamOption::BINARY);

        if (useCoherent)
        {
            #ifdef USE_COHERENT
            streamOpt.format(IOstreamOption::COHERENT);
            runTime.writeFormat(IOstreamOption::COHERENT);

            mesh.writeObject(streamOpt, true);

            Info<< nl
                << "Specified -coherent (instance: "
                << mesh.pointsInstance() << ")" << endl;

            const auto& coherent = CoherentMesh::New(mesh);

            Info<< "    points  = "
                << coherent.globalPointOffsets().totalSize() << nl
                << "    cells   = "
                << coherent.globalCellOffsets().totalSize() << nl
                << "    patches = "
                << coherent.nNonProcessorPatches() << nl;

            #else
            Info<< "Warning: -coherent ignored" << nl;
            #endif
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

            for (const auto& fld : fields)
            {
                fld.regIOobject::writeObject(streamOpt, writeOnProc);
            }
        }

        if (useCoherent)
        {
            #ifdef USE_COHERENT
            SliceStreamRepo::closeInstance();
            #endif
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
    }
    else
    {
        label fieldSize = 1000;

        if (args.readIfPresent("global", fieldSize))
        {
            fieldSize /= nProcsEff;
        }
        else
        {
            args.readIfPresent("local", fieldSize);
        }

        if (!writeOnProc)
        {
            fieldSize = 0;
        }

        const globalIndex giCells(fieldSize);

        // Create fields
        Info<< nl << "Create " << nFields << " fields" << nl
            << "field-size:" << fieldSize
            << " total-size:" << giCells.totalSize() << nl;


        PtrList<IOField<scalar>> fields(nFields);

        {
            IOobject io
            (
                "field",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            forAll(fields, fieldi)
            {
                io.resetHeader("field" + Foam::name(fieldi));

                fields.set
                (
                    fieldi,
                    new IOField<scalar>(io, fieldSize)
                );

                auto& fld = fields[fieldi];
                std::iota(fld.begin(), fld.end(), scalar(0));
            }
        }

        IOstreamOption streamOpt(IOstreamOption::BINARY);

        if (useCoherent)
        {
            Info<< "Warning: -coherent ignored" << nl;
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

            for (const auto& fld : fields)
            {
                fld.regIOobject::writeObject(streamOpt, writeOnProc);
            }
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
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
