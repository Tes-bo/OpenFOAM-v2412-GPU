/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "commSchedule.H"
#include "IOstreams.H"
#include "IOmanip.H"
#include "StringStream.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(commSchedule, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Count the number of outstanding communications for a single processor
static label outstandingComms
(
    const labelUList& commToSchedule,
    const DynamicList<label>& procComms
)
{
    label nOutstanding = 0;

    for (const label commPairi : procComms)
    {
        if (commToSchedule[commPairi] == -1)
        {
            ++nOutstanding;
        }
    }
    return nOutstanding;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::commSchedule::commSchedule
(
    const label nProcs,
    const List<labelPair>& comms
)
:
    schedule_(comms.size()),
    procSchedule_(nProcs)
{
    // Determine comms per processor.
    List<DynamicList<label>> procToComms(nProcs);

    forAll(comms, commPairi)
    {
        const label proc0 = comms[commPairi].first();
        const label proc1 = comms[commPairi].second();

        if (proc0 < 0 || proc0 >= nProcs || proc1 < 0 || proc1 >= nProcs)
        {
            FatalErrorInFunction
                << "Illegal processor(s): "
                << comms[commPairi] << abort(FatalError);
        }

        procToComms[proc0].push_back(commPairi);
        procToComms[proc1].push_back(commPairi);
    }
    // Note: no need to shrink procToComms. Are small.

    if (debug && UPstream::master())
    {
        Pout<< "commSchedule : Wanted communication:" << endl;

        forAll(comms, i)
        {
            const labelPair& twoProcs = comms[i];

            Pout<< i << ": "
                << twoProcs.first() << " <-> " << twoProcs.second() << endl;
        }
        Pout<< endl;


        Pout<< "commSchedule : Schedule:" << endl;

        // Print header. Use buffered output to prevent parallel output messing
        // up.
        {
            OStringStream os;
            os  << "iter|";
            for (int i = 0; i < nProcs; i++)
            {
                os  << setw(3) << i;
            }
            Pout<< os.str().c_str() << endl;
        }
        {
            OStringStream os;
            os  << "----+";
            for (int i = 0; i < nProcs; i++)
            {
                os  << "---";
            }
            Pout<< os.str().c_str() << endl;
        }
    }

    // Schedule all. Note: crap scheduler. Assumes all communication takes
    // equally long.

    label nScheduled = 0;

    label iter = 0;

    // Per index into comms the time when it was scheduled
    labelList commToSchedule(comms.size(), -1);

    while (nScheduled < comms.size())
    {
        label oldNScheduled = nScheduled;

        // Find unscheduled comms. This is the comms where the two processors
        // still have the most unscheduled comms.

        boolList busy(nProcs, false);

        while (true)
        {
            label maxComm = -1;
            label maxNeed = labelMin;

            forAll(comms, commPairi)
            {
                const label proc0 = comms[commPairi].first();
                const label proc1 = comms[commPairi].second();

                if
                (
                    commToSchedule[commPairi] == -1  // unscheduled
                 && !busy[proc0]
                 && !busy[proc1]
                )
                {
                    label need =
                    (
                        outstandingComms(commToSchedule, procToComms[proc0])
                      + outstandingComms(commToSchedule, procToComms[proc1])
                    );

                    if (maxNeed < need)
                    {
                        maxNeed = need;
                        maxComm = commPairi;
                    }
                }
            }


            if (maxComm == -1)
            {
                // Found no unscheduled procs.
                break;
            }

            // Schedule commPairi in this iteration
            commToSchedule[maxComm] = nScheduled++;
            busy[comms[maxComm].first()] = true;
            busy[comms[maxComm].second()] = true;
        }

        if (debug && UPstream::master())
        {
            label nIterComms = nScheduled-oldNScheduled;

            if (nIterComms > 0)
            {
                labelList procToComm(nProcs, -1);

                forAll(commToSchedule, commPairi)
                {
                    const label sched = commToSchedule[commPairi];

                    if (sched >= oldNScheduled && sched < nScheduled)
                    {
                        const label proc0 = comms[commPairi].first();
                        const label proc1 = comms[commPairi].second();
                        procToComm[proc0] = commPairi;
                        procToComm[proc1] = commPairi;
                    }
                }

                // Print it
                OStringStream os;
                os  << setw(3) << iter << " |";
                forAll(procToComm, proci)
                {
                    if (procToComm[proci] == -1)
                    {
                        os  << "   ";
                    }
                    else
                    {
                        os  << setw(3) << procToComm[proci];
                    }
                }
                Pout<< os.str().c_str() << endl;
            }
        }

        iter++;
    }

    if (debug && UPstream::master())
    {
        Pout<< endl;
    }


    // Sort commToSchedule to obtain order in comms

    Foam::sortedOrder(commToSchedule, schedule_);

    // Sort schedule_ by processor

    labelList nProcScheduled(nProcs, Zero);

    // Count
    for (const label commPairi : schedule_)
    {
        const labelPair& twoProcs = comms[commPairi];

        nProcScheduled[twoProcs.first()]++;
        nProcScheduled[twoProcs.second()]++;
    }

    // Allocate
    forAll(procSchedule_, proci)
    {
        procSchedule_[proci].resize_nocopy(nProcScheduled[proci]);
    }

    nProcScheduled = 0;

    // Fill
    for (const label commPairi : schedule_)
    {
        const labelPair& twoProcs = comms[commPairi];

        const label proc0 = twoProcs.first();
        const label proc1 = twoProcs.second();

        procSchedule_[proc0][nProcScheduled[proc0]++] = commPairi;
        procSchedule_[proc1][nProcScheduled[proc1]++] = commPairi;
    }

    if (debug && UPstream::master())
    {
        Pout<< "commSchedule::commSchedule : Per processor:" << endl;

        forAll(procSchedule_, proci)
        {
            const labelList& procComms = procSchedule_[proci];

            Pout<< "Processor " << proci << " talks to processors:" << endl;

            for (const label commPairi : procComms)
            {
                const labelPair& twoProcs = comms[commPairi];

                Pout<< "    "
                    << (proci == twoProcs[1] ? twoProcs[0] : twoProcs[1])
                    << endl;
            }
        }
        Pout<< endl;
    }
}


// ************************************************************************* //
