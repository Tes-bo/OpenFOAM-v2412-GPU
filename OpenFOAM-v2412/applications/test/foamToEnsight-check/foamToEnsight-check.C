/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    foamToEnsight-check

Description
    Check data sizes for conversion to ensight format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "Pstream.H"
#include "HashOps.H"
#include "regionProperties.H"

#include "fvc.H"
#include "faMesh.H"
#include "fvMesh.H"

// file-format/conversion
#include "ensightFaMesh.H"
#include "ensightMesh.H"

using namespace Foam;


void printStats(const FixedList<label, 3>& stats, const char *what = "")
{
    Info<< what << "max-comm: "<< stats[0] << nl
        << what << "max-size: "<< stats[1] << nl
        << what << "off-proc: "<< stats[2] << nl;
}


template<class EnsightPartType>
FixedList<label, 3> printPartInfo
(
    const EnsightPartType& part,
    int verbose = 0
)
{
    Info<< "part: " << part.name().c_str() << nl
        << "    size: "
        << (Pstream::parRun() ? part.total() : part.size())
        << " (";

    FixedList<label, 3> stats(Zero);

    label& maxComm = stats[0];
    label& maxSize = stats[1];
    label& totNonLocalSize = stats[2];

    for (int typei=0; typei < EnsightPartType::nTypes; ++typei)
    {
        const auto etype = typename EnsightPartType::elemType(typei);

        if (typei) Info<< ' ';
        Info<< EnsightPartType::elemNames[etype] << ": "
            << (Pstream::parRun() ? part.total(etype) : part.size(etype));

        label elemCount = part.size(etype);
        label commCount = (Pstream::master() ? label(0) : elemCount);
        label nonLocalCount = commCount;

        if (Pstream::parRun())
        {
            reduce(elemCount, maxOp<label>());
            reduce(commCount, maxOp<label>());

            reduce(nonLocalCount, sumOp<label>());
        }

        maxComm = max(maxComm, commCount);
        maxSize = max(maxSize, elemCount);
        totNonLocalSize = max(totNonLocalSize, nonLocalCount);
    }
    Info<< ")" << endl;

    if (verbose && Pstream::parRun() && part.total())
    {
        for (int typei=0; typei < EnsightPartType::nTypes; ++typei)
        {
            const auto etype = typename EnsightPartType::elemType(typei);

            label elemCount = part.size(etype);
            label totCount = part.total(etype);

            Info<< "    "
                << EnsightPartType::elemNames[etype] << ": "
                << totCount;

            if (totCount)
            {
                labelList sizes(UPstream::listGatherValues(elemCount));

                Info<< "  ";
                sizes.writeList(Info);
            }

            Info<< endl;
        }
    }

    printStats(stats, "    ");

    return stats;
}


void printInfo(const ensightMesh& mesh, int verbose = 0)
{
    FixedList<label, 3> cellStats(Zero);
    FixedList<label, 3> faceStats(Zero);

    for (const auto& iter : mesh.cellZoneParts().csorted())
    {
        FixedList<label, 3> stats = printPartInfo(iter.val(), verbose);

        for (label i=0; i < 3; ++i)
        {
            cellStats[i] = max(cellStats[i], stats[i]);
        }
    }

    for (const auto& iter : mesh.faceZoneParts().csorted())
    {
        FixedList<label, 3> stats = printPartInfo(iter.val(), verbose);

        for (label i=0; i < 3; ++i)
        {
            faceStats[i] = max(faceStats[i], stats[i]);
        }
    }

    for (const auto& iter : mesh.boundaryParts().csorted())
    {
        FixedList<label, 3> stats = printPartInfo(iter.val(), verbose);

        for (label i=0; i < 3; ++i)
        {
            faceStats[i] = max(faceStats[i], stats[i]);
        }
    }

    Info<< nl
        << "===============" << nl;
    printStats(cellStats, "cell ");

    Info<< nl;
    printStats(faceStats, "face ");

    Info<< "===============" << endl;
}


void printInfo(const ensightFaMesh& mesh, int verbose = 0)
{
    printPartInfo(mesh.areaPart());
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Check data sizes for conversion of OpenFOAM to Ensight format"
    );
    // timeSelector::addOptions();

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");

    argList::addVerboseOption();

    #include "addAllRegionOptions.H"

    argList::addBoolOption
    (
        "no-boundary",  // noPatches
        "Suppress writing any patches"
    );
    argList::addBoolOption
    (
        "no-internal",
        "Suppress writing the internal mesh"
    );
    argList::addBoolOption
    (
        "no-cellZones",
        "Suppress writing any cellZones"
    );
    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress output of finite-area mesh/fields",
        true  // mark as an advanced option
    );

    #include "setRootCase.H"

    // ------------------------------------------------------------------------
    // Configuration

    const int optVerbose = args.verbose();
    const bool doBoundary    = !args.found("no-boundary");
    const bool doInternal    = !args.found("no-internal");
    const bool doCellZones   = !args.found("no-cellZones");
    const bool doFiniteArea  = !args.found("no-finite-area");

    ensightMesh::options writeOpts;
    writeOpts.useBoundaryMesh(doBoundary);
    writeOpts.useInternalMesh(doInternal);
    writeOpts.useCellZones(doCellZones);

    // ------------------------------------------------------------------------

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    // Handle -allRegions, -regions, -region
    #include "getAllRegionOptions.H"

    // ------------------------------------------------------------------------

    #include "createNamedMeshes.H"

    // ------------------------------------------------------------------------
    /// #include "createMeshAccounting.H"

    PtrList<ensightMesh> ensightMeshes(regionNames.size());
    PtrList<faMesh> meshesFa(regionNames.size());
    PtrList<ensightFaMesh> ensightMeshesFa(regionNames.size());

    forAll(regionNames, regioni)
    {
        const fvMesh& mesh = meshes[regioni];

        ensightMeshes.set
        (
            regioni,
            new ensightMesh(mesh, writeOpts)
        );
        ensightMeshes[regioni].verbose(optVerbose);


        if (doFiniteArea)
        {
            autoPtr<faMesh> faMeshPtr(faMesh::TryNew(mesh));

            if (faMeshPtr)
            {
                meshesFa.set(regioni, std::move(faMeshPtr));

                ensightMeshesFa.set
                (
                    regioni,
                    new ensightFaMesh(meshesFa[regioni])
                );
                ensightMeshesFa[regioni].verbose(optVerbose);
            }
        }
    }

    // ------------------------------------------------------------------------


    if (Pstream::master())
    {
        Info<< "Checking " << timeDirs.size() << " time steps" << nl;
    }

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            // const word& regionDir = polyMesh::regionName(regionName);

            auto& mesh = meshes[regioni];

            polyMesh::readUpdateState meshState = mesh.readUpdate();
            const bool moving = (meshState != polyMesh::UNCHANGED);

            auto& ensMesh = ensightMeshes[regioni];

            // Finite-area (can be missing)
            auto* ensFaMeshPtr = ensightMeshesFa.get(regioni);

            if (moving)
            {
                ensMesh.expire();
                ensMesh.correct();

                if (ensFaMeshPtr)
                {
                    ensFaMeshPtr->expire();
                    ensFaMeshPtr->correct();
                }
            }

            if (moving || timei == 0)  // report
            {
                if (regionNames.size() > 1)
                {
                    Info<< "region=" << regionName << nl;
                }

                printInfo(ensMesh, optVerbose);

                if (ensFaMeshPtr)
                {
                    printInfo(*ensFaMeshPtr, optVerbose);
                }
            }
        }
    }

    Info<< "\nEnd"<< nl << endl;

    return 0;
}


// ************************************************************************* //
