/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "multiWorldConnectionsObject.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiWorldConnections, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Combine world-to-world connections.
// Forward connection = 1, Backward connection = 2, Both = 3
struct worldConnectBitOrEq
{
    void operator()(EdgeMap<unsigned>& a, const EdgeMap<unsigned>& b) const
    {
        forAllConstIters(b, iter)
        {
            a(iter.key()) |= iter.val();
        }
    }
};


static void printDOT(Ostream& os, const EdgeMap<unsigned>& connections)
{
    os << nl << "// Multiworld communication graph:" << nl;
    os.beginBlock("graph");

    // Graph Nodes == worlds
    label worldi = 0;
    for (const word& worldName : UPstream::allWorlds())
    {
        os.indent();
        os  << worldi << " [xlabel=" << worldi
            << ",label=\"" << worldName << "\"]" << nl;

        ++worldi;
    }
    os  << nl;

    // Graph Edges == connections
    for (const edge& connect : connections.sortedToc())
    {
        os.indent();
        os  << connect.first() << " -- " << connect.second();

        // Mismatched forward/backward connections?
        if (connections.lookup(connect, 0u) != 3u)
        {
            os << " [style=dashed] // mismatched?";
        }
        os << nl;
    }

    os.endBlock();

    os << "// end graph" << nl;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::edge Foam::multiWorldConnections::worldPair(const label otherWorld)
{
    if (otherWorld < 0 || !Pstream::parRun())
    {
        Perr<< "ignore: no world or non-parallel" << endl;
        return edge(-1, -1);
    }
    else if (UPstream::allWorlds().size() <= otherWorld)
    {
        Perr<< "ignore: invalid world: " << otherWorld << endl;
        return edge(-1, -1);
    }

    const label thisWorldID = UPstream::myWorldID();

    // The worlds (sorted)
    return edge(thisWorldID, otherWorld, true);
}


Foam::edge Foam::multiWorldConnections::worldPair(const word& otherWorld)
{
    if (otherWorld.empty() || !Pstream::parRun())
    {
        Perr<< "ignore: no world or non-parallel" << endl;
        return edge(-1, -1);
    }

    const label thisWorldID = UPstream::myWorldID();
    const label otherWorldID = UPstream::allWorlds().find(otherWorld);

    if (otherWorldID < 0)
    {
        FatalErrorInFunction
            << "Cannot find world " << otherWorld
            << " in set of worlds " << flatOutput(UPstream::allWorlds())
            << exit(FatalError);
    }

    // The worlds (sorted)
    return edge(thisWorldID, otherWorldID, true);
}


Foam::label Foam::multiWorldConnections::createCommunicator(const edge& worlds)
{
    // Fallback: do not create, just use local world
    label comm = UPstream::worldComm;

    if (!worlds.good())
    {
        return comm;
    }

    const labelList& worldIDs = UPstream::worldIDs();

    DynamicList<label> subRanks(worldIDs.size());
    forAll(worldIDs, proci)
    {
        if (worlds.found(worldIDs[proci]))
        {
            subRanks.push_back(proci);
        }
    }

    // Allocate new communicator with global world
    comm = UPstream::allocateCommunicator(UPstream::commGlobal(), subRanks);

    if (debug & 2)
    {
        Pout<< "multiWorld::communicator :"
            << " between " << UPstream::allWorlds()[worlds.first()]
            << " and " << UPstream::allWorlds()[worlds.second()]
            << " sub-ranks: " << subRanks
            << " comm:" << comm << endl;
    }

    return comm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiWorldConnections::multiWorldConnections(const Time& runTime)
:
    MeshObject_type(runTime)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::multiWorldConnections&
Foam::multiWorldConnections::New(const Time& runTime)
{
    return MeshObject_type::New(runTime);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiWorldConnections::empty() const noexcept
{
    return table_.empty();
}


Foam::label Foam::multiWorldConnections::size() const noexcept
{
    return table_.size();
}


void Foam::multiWorldConnections::createComms()
{
    // Need new communicator(s)

    const label thisWorldID = UPstream::myWorldID();

    EdgeMap<unsigned> allConnections;
    forAllConstIters(table_, iter)
    {
        const edge& connect = iter.key();

        allConnections.insert
        (
            connect,
            (connect.first() == thisWorldID ? 1u : 2u)
        );
    }


    // Use MPI_COMM_WORLD
    const label oldWarnComm = UPstream::commWarn(UPstream::commGlobal());
    const label oldWorldComm = UPstream::commWorld(UPstream::commGlobal());

    Pstream::combineReduce(allConnections, worldConnectBitOrEq());

    // Check for mismatched connections
    label brokenConnections = 0;

    forAllConstIters(allConnections, iter)
    {
        // Mismatched forward/backward connections?
        if (iter.val() != 3u)
        {
            ++brokenConnections;
        }
    }

    if (brokenConnections)
    {
        // Restore communicator settings
        UPstream::commWarn(oldWarnComm);
        UPstream::commWorld(oldWorldComm);

        FatalErrorInFunction
            << "Has " << brokenConnections
            << " broken world-world connections";

        printDOT(FatalError, allConnections);

        FatalError << exit(FatalError);
    }
    else
    {
        // NOTE: process in sorted order to ensure proper
        // synchronization on all worlds and all processors

        for (const edge& connect : allConnections.sortedToc())
        {
            // Process known connections without communicators.
            // - create a communicator and cache its value

            auto iter = table_.find(connect);
            if (iter.good() && iter.val() == -1)
            {
                iter.val() = createCommunicator(connect);
            }
        }

        // Restore communicator settings
        UPstream::commWarn(oldWarnComm);
        UPstream::commWorld(oldWorldComm);
    }

    if (debug)
    {
        printDOT(Info, allConnections);
    }
}


bool Foam::multiWorldConnections::addConnectionById(const label otherWorld)
{
    // The worlds (sorted)
    edge worlds(worldPair(otherWorld));

    if (!worlds.good())
    {
        return false;
    }

    const bool added = table_.insert(worlds, -1);

    Pout<< (added ? "Add" : "Existing") << " connection from "
        << UPstream::myWorld() << " to " << otherWorld << nl;

    return added;
}


bool Foam::multiWorldConnections::addConnectionByName(const word& otherWorld)
{
    // The worlds (sorted)
    edge worlds(worldPair(otherWorld));

    if (!worlds.good())
    {
        return false;
    }

    const bool added = table_.insert(worlds, -1);

    Pout<< (added ? "Add" : "Existing") << " connection from "
        << UPstream::myWorld() << " to " << otherWorld << nl;

    return added;
}


Foam::label Foam::multiWorldConnections::getCommById
(
    const label otherWorldID
) const
{
    // Default: use local world
    label comm = UPstream::worldComm;

    // The communication worlds (sorted)
    edge worlds(worldPair(otherWorldID));

    if (!worlds.good())
    {
        return comm;
    }

    const auto iter = table_.cfind(worlds);

    if (!iter.good())
    {
        FatalErrorInFunction
            << "No connection registered for worlds " << worlds
            << exit(FatalError);
    }

    // Get cached value, or allocate ALL known communicators
    comm = iter.val();

    if (comm == -1)
    {
        // Need new communicator(s)
        const_cast<multiWorldConnections&>(*this).createComms();

        // Retrieve from table cache
        comm = table_.lookup(worlds, UPstream::worldComm);
    }

    return comm;
}


Foam::label Foam::multiWorldConnections::getCommByName
(
    const word& otherWorld
) const
{
    // Default: use local world
    label comm = UPstream::worldComm;

    // The communication worlds (sorted)
    edge worlds(worldPair(otherWorld));

    if (!worlds.good())
    {
        return comm;
    }

    const auto iter = table_.cfind(worlds);

    if (!iter.good())
    {
        FatalErrorInFunction
            << "No connection registered for worlds " << worlds
            << exit(FatalError);
    }

    // Get cached value, or allocate ALL known communicators
    comm = iter.val();

    if (comm == -1)
    {
        // Need new communicator(s)
        const_cast<multiWorldConnections&>(*this).createComms();

        // Retrieve from table cache
        comm = table_.lookup(worlds, UPstream::worldComm);
    }

    return comm;
}


Foam::labelList Foam::multiWorldConnections::comms() const
{
    labelList list(table_.size());

    if (list.empty())
    {
        // Default: use local world
        list.resize(1, UPstream::worldComm);
    }
    else
    {
        forAllConstIters(table_, iter)
        {
            if (iter.val() == -1)
            {
                // Need new communicator(s)
                const_cast<multiWorldConnections&>(*this).createComms();
                break;
            }
        }

        // Retrieve values from table cache
        label i = 0;

        forAllConstIters(table_, iter)
        {
            list[i] = iter.val();
            ++i;
        }

        Foam::sort(list);  // Consistent order!
    }

    return list;
}


// ************************************************************************* //
