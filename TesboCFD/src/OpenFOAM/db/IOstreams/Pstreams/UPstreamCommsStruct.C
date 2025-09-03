/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "UPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::commsStruct::commsStruct
(
    const label above,
    labelList&& below,
    labelList&& allBelow,
    labelList&& allNotBelow
)
:
    above_(above),
    below_(std::move(below)),
    allBelow_(std::move(allBelow)),
    allNotBelow_(std::move(allNotBelow))
{}


Foam::UPstream::commsStruct::commsStruct
(
    const label numProcs,
    const label myProcID,
    const label above,
    const labelUList& below,
    const labelUList& allBelow
)
:
    above_(above),
    below_(below),
    allBelow_(allBelow),
    allNotBelow_(numProcs - allBelow.size() - 1)
{
    List<bool> isNotBelow(numProcs, true);

    // Exclude self
    isNotBelow[myProcID] = false;

    // Exclude allBelow
    for (const label proci : allBelow)
    {
        isNotBelow[proci] = false;
    }

    // Compacting to obtain allNotBelow_
    label nNotBelow = 0;
    forAll(isNotBelow, proci)
    {
        if (isNotBelow[proci])
        {
            allNotBelow_[nNotBelow++] = proci;
        }
    }

    if (nNotBelow != allNotBelow_.size())
    {
        FatalErrorInFunction
            << "Problem: " << nNotBelow << " != " << allNotBelow_.size() << nl
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// This outputs as depth-first, but graphviz sorts that for us
void Foam::UPstream::commsStruct::printGraph
(
    Ostream& os,
    const UList<UPstream::commsStruct>& comms,
    const label proci
)
{
    // if (proci >= comms.size()) return;  // Extreme safety!

    const auto& below = comms[proci].below();

    if (proci == 0)
    {
        os << nl << "// communication graph:" << nl;
        os.beginBlock("graph");

        if (below.empty())
        {
            // A graph with a single-node (eg, self-comm)
            os << indent << proci << nl;
        }
    }

    int pos = 0;

    for (const label nbrProci : below)
    {
        if (pos)
        {
            os << "  ";
        }
        else
        {
            os << indent;
        }
        os << proci << " -- " << nbrProci;

        if (++pos >= 4)  // Max 4 items per line
        {
            pos = 0;
            os << nl;
        }
    }

    if (pos)
    {
        os << nl;
    }

    for (const label nbrProci : below)
    {
        // if (proci == nbrProci) continue;  // Extreme safety!
        printGraph(os, comms, nbrProci);
    }

    if (proci == 0)
    {
        os.endBlock();

        os << "// end graph" << nl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UPstream::commsStruct::nProcs() const
{
    return (1 + allBelow_.size() + allNotBelow_.size());
}


void Foam::UPstream::commsStruct::reset()
{
    above_ = -1;
    below_.clear();
    allBelow_.clear();
    allNotBelow_.clear();
}


void Foam::UPstream::commsStruct::reset
(
    const label procID,
    const label numProcs
)
{
    reset();

    label above(-1);
    DynamicList<label> below;
    DynamicList<label> allBelow;

    if (numProcs < UPstream::nProcsSimpleSum)
    {
        // Linear schedule

        if (procID == 0)
        {
            below = identity(numProcs-1, 1);
            allBelow = below;
        }
        else
        {
            above = 0;
        }
    }
    else
    {
        // Use tree like schedule. For 8 procs:
        // (level 0)
        //      0 receives from 1
        //      2 receives from 3
        //      4 receives from 5
        //      6 receives from 7
        // (level 1)
        //      0 receives from 2
        //      4 receives from 6
        // (level 2)
        //      0 receives from 4
        //
        // The sends/receives for all levels are collected per processor
        // (one send per processor; multiple receives possible) creating
        // a table:
        //
        // So per processor:
        // proc     receives from   sends to
        // ----     -------------   --------
        //  0       1,2,4           -
        //  1       -               0
        //  2       3               0
        //  3       -               2
        //  4       5               0
        //  5       -               4
        //  6       7               4
        //  7       -               6

        label mod = 0;

        for (label step = 1; step < numProcs; step = mod)
        {
            mod = step * 2;

            if (procID % mod)
            {
                above = procID - (procID % mod);
                break;
            }
            else
            {
                for
                (
                    label j = procID + step;
                    j < numProcs && j < procID + mod;
                    j += step
                )
                {
                    below.push_back(j);
                }
                for
                (
                    label j = procID + step;
                    j < numProcs && j < procID + mod;
                    j++
                )
                {
                    allBelow.push_back(j);
                }
            }
        }
    }

    *this = UPstream::commsStruct(numProcs, procID, above, below, allBelow);
}


// * * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * //

template<>
Foam::UPstream::commsStruct&
Foam::UList<Foam::UPstream::commsStruct>::operator[](const label procID)
{
    auto& val = this->v_[procID];   // or this->data()[procID]

    if (val.nProcs() != size())
    {
        // Create/update
        val.reset(procID, size());
    }

    return val;
}


template<>
const Foam::UPstream::commsStruct&
Foam::UList<Foam::UPstream::commsStruct>::operator[](const label procID) const
{
    return const_cast<UList<UPstream::commsStruct>&>(*this).operator[](procID);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::UPstream::commsStruct::operator==(const commsStruct& comm) const
{
    return
    (
        (above_ == comm.above())
     && (below_ == comm.below())
     // && (allBelow_ == comm.allBelow())
     // && (allNotBelow_ == comm.allNotBelow())
    );
}


bool Foam::UPstream::commsStruct::operator!=(const commsStruct& comm) const
{
    return !operator==(comm);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const UPstream::commsStruct& comm)
{
    os  << comm.above() << nl << token::SPACE << token::SPACE;
    comm.below().writeList(os) << nl << token::SPACE << token::SPACE;
    comm.allBelow().writeList(os) << nl << token::SPACE << token::SPACE;
    comm.allNotBelow().writeList(os);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
