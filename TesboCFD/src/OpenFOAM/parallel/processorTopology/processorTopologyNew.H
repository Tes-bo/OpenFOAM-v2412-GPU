/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

InClass
    Foam::processorTopology

Description
    Define the processor-processor connection table by walking a list of
    patches and detecting the processor ones.

    This has been split off as a separate include from processorTopology
    to minimize dependencies.

Warning
    Does not currently correctly support multiple processor
    patches connecting two processors.

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorTopologyNew_H
#define Foam_processorTopologyNew_H

#include "processorTopology.H"
#include "commSchedule.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class ProcPatch, class PatchListType>
Foam::processorTopology Foam::processorTopology::New
(
    const PatchListType& patches,
    const label comm
)
{
    processorTopology topo;

    topo.comm_ = comm;  // The communicator
    auto& schedule = topo.patchSchedule_;

    schedule.resize(2*patches.size());


    // The evaluation number within the schedule
    label patchEvali = 0;

    // Number of processor patches
    label numProcPatches = 0;

    //
    // 1. Schedule all non-processor patches
    //

    forAll(patches, patchi)
    {
        if (isA<ProcPatch>(patches[patchi]))
        {
            ++numProcPatches;
        }
        else
        {
            schedule[patchEvali++].setInitEvaluate(patchi);
            schedule[patchEvali++].setEvaluate(patchi);
        }
    }


    // Assemble processor patch information
    if (UPstream::parRun() && numProcPatches)
    {
        // Create reverse map (from proc to patch)
        // - assumes single connections between processors

        auto& patchMap = topo.procPatchMap_;
        patchMap.reserve(numProcPatches);

        forAll(patches, patchi)
        {
            const auto* cpp = isA<ProcPatch>(patches[patchi]);
            if (cpp)
            {
                patchMap.set(cpp->neighbProcNo(), patchi);
            }
        }
    }


    //
    // 2. Handle processor patches
    //

    if
    (
        UPstream::parRun()
     && UPstream::defaultCommsType == UPstream::commsTypes::scheduled
    )
    {
        const label myProci = UPstream::myProcNo(comm);
        const label nProcs = UPstream::nProcs(comm);

        // Synchronized on all processors
        const auto& procToProcs = topo.procAdjacency();

        // Determine the schedule for all processor patches.
        // Insert processor pair once to determine the schedule.
        // Each processor pair stands for both send and receive.

        label nComms = 0;
        for (const labelList& nbrProcs : procToProcs)
        {
            nComms += nbrProcs.size();
        }
        DynamicList<labelPair> comms(nComms/2);

        forAll(procToProcs, proci)
        {
            for (const label nbrProci : procToProcs[proci])
            {
                if (proci < nbrProci)
                {
                    // Owner to neighbour connection
                    comms.push_back(labelPair(proci, nbrProci));
                }
            }
        }

        // Determine a schedule.

        labelList mySchedule
        (
            commSchedule(nProcs, comms).procSchedule()[myProci]
        );

        for (const label scheduleIndex : mySchedule)
        {
            // Get the other processor
            label nbrProci = comms[scheduleIndex].first();
            if (nbrProci == myProci)
            {
                nbrProci = comms[scheduleIndex].second();
            }
            const label patchi = topo.procPatchLookup(nbrProci);

            if (myProci > nbrProci)
            {
                schedule[patchEvali++].setInitEvaluate(patchi);
                schedule[patchEvali++].setEvaluate(patchi);
            }
            else
            {
                schedule[patchEvali++].setEvaluate(patchi);
                schedule[patchEvali++].setInitEvaluate(patchi);
            }
        }
    }
    else
    {
        // Non-blocking schedule for processor patches

        if (numProcPatches)
        {
            forAll(patches, patchi)
            {
                if (isA<ProcPatch>(patches[patchi]))
                {
                    schedule[patchEvali].setInitEvaluate(patchi);
                    schedule[patchEvali + numProcPatches].setEvaluate(patchi);
                    ++patchEvali;
                }
            }
        }
    }

    return topo;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
