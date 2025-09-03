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

\*---------------------------------------------------------------------------*/

#include "processorTopology.H"
#include "Pstream.H"   // For allGatherList()

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorTopology::processorTopology() noexcept
:
    comm_(UPstream::worldComm)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::processorTopology::procNeighbours() const
{
    if (procNeighbours_.empty() && !procPatchMap_.empty())
    {
        // My neighbouring procs in ascending sorted order
        procNeighbours_ = procPatchMap_.sortedToc();
    }

    return procNeighbours_;
}


// May be useful in the future...
// ------------------------------
//
// const Foam::labelUList Foam::processorTopology::below() const
// {
//     const auto& all = procNeighbours();
//
//     const auto* pivot = std::upper_bound
//     (
//         all.begin(),
//         all.end(),
//         UPstream::myProcNo(comm_)
//     );
//
//     if (pivot != all.end())
//     {
//         return UList<label>
//         (
//             const_cast<label*>(all.begin()),
//             (pivot - all.begin())
//         );
//     }
//     return UList<label>();
// }
//
//
// const Foam::labelUList Foam::processorTopology::above() const
// {
//     const auto& all = procNeighbours();
//
//     const auto* pivot = std::upper_bound
//     (
//         all.begin(),
//         all.end(),
//         UPstream::myProcNo(comm_)
//     );
//     if (pivot != all.end())
//     {
//         return UList<label>
//         (
//             const_cast<label*>(pivot),
//             (all.end() - pivot)
//         );
//     }
//     return UList<label>();
// }


const Foam::labelListList& Foam::processorTopology::procAdjacency() const
{
    if (UPstream::parRun() && procAdjacencyTable_.empty())
    {
        procAdjacencyTable_.resize(UPstream::nProcs(comm_));

        // My neighbouring procs in ascending sorted order
        procAdjacencyTable_[UPstream::myProcNo(comm_)]
            = procPatchMap_.sortedToc();

        // Synchronize on all processors
        Pstream::allGatherList(procAdjacencyTable_, UPstream::msgType(), comm_);
    }

    return procAdjacencyTable_;
}


// ************************************************************************* //
