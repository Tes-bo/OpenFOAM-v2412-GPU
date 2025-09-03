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

\*---------------------------------------------------------------------------*/

#include "hostUncollatedFileOperation.H"
#include "fileOperationInitialise.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(hostUncollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostUncollatedFileOperation,
        word
    );
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostUncollatedFileOperation,
        comm
    );

    // Threaded MPI: not required
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        fileOperationInitialise_unthreaded,
        word,
        hostUncollated
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

void Foam::fileOperations::hostUncollatedFileOperation::init(bool verbose)
{
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (verbose)
    {
        DetailInfo
            << "I/O    : " << this->type() << nl;

        if (ioRanks_.size())
        {
            fileOperation::printRanks();
        }
    }
}


Foam::fileOperations::hostUncollatedFileOperation::hostUncollatedFileOperation
(
    bool verbose
)
:
    masterUncollatedFileOperation
    (
        getCommPattern(),
        false,  // distributedRoots
        false   // verbose
    ),
    managedComm_(getManagedComm(comm_))  // Possibly locally allocated
{
    init(verbose);
}


Foam::fileOperations::hostUncollatedFileOperation::hostUncollatedFileOperation
(
    const Tuple2<label, labelList>& commAndIORanks,
    const bool distributedRoots,
    bool verbose
)
:
    masterUncollatedFileOperation
    (
        commAndIORanks,
        distributedRoots,
        false   // verbose
    ),
    managedComm_(-1)  // Externally managed
{
    init(verbose);
}


void Foam::fileOperations::hostUncollatedFileOperation::storeComm() const
{
    // From externally -> locally managed
    managedComm_ = getManagedComm(comm_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::hostUncollatedFileOperation::
~hostUncollatedFileOperation()
{
    // Wait for any outstanding file operations
    flush();

    UPstream::freeCommunicator(managedComm_);
}


// ************************************************************************* //
