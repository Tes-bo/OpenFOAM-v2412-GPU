/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "hostCollatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(hostCollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostCollatedFileOperation,
        word
    );
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostCollatedFileOperation,
        comm
    );

    // Threaded MPI: depending on buffering
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        fileOperationInitialise_collated,
        word,
        hostCollated
    );
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
// Construction helper: self/world/local communicator and IO ranks
static Tuple2<label, labelList> getCommPattern()
{
    // Default is COMM_WORLD (single master)
    Tuple2<label, labelList> commAndIORanks
    (
        UPstream::worldComm,
        fileOperation::getGlobalIORanks()
    );

    if (commAndIORanks.second().empty())
    {
        // Default: one master per host
        commAndIORanks.second() = fileOperation::getGlobalHostIORanks();
    }

    if (UPstream::parRun() && commAndIORanks.second().size() > 1)
    {
        // Multiple masters: ranks for my IO range
        commAndIORanks.first() = UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            fileOperation::subRanks(commAndIORanks.second())
        );
    }

    return commAndIORanks;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::fileOperations::hostCollatedFileOperation::init(bool verbose)
{
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (verbose)
    {
        this->printBanner(ioRanks_.size());
    }
}


Foam::fileOperations::hostCollatedFileOperation::hostCollatedFileOperation
(
    bool verbose
)
:
    collatedFileOperation
    (
        getCommPattern(),
        false,  // distributedRoots
        false   // verbose
    ),
    managedComm_(getManagedComm(comm_))  // Possibly locally allocated
{
    init(verbose);
}


Foam::fileOperations::hostCollatedFileOperation::hostCollatedFileOperation
(
    const Tuple2<label, labelList>& commAndIORanks,
    const bool distributedRoots,
    bool verbose
)
:
    collatedFileOperation
    (
        commAndIORanks,
        distributedRoots,
        false   // verbose
    ),
    managedComm_(-1)  // Externally managed
{
    init(verbose);
}


void Foam::fileOperations::hostCollatedFileOperation::storeComm() const
{
    // From externally -> locally managed
    managedComm_ = getManagedComm(comm_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::hostCollatedFileOperation::~hostCollatedFileOperation()
{
    UPstream::freeCommunicator(managedComm_);
}


// ************************************************************************* //
