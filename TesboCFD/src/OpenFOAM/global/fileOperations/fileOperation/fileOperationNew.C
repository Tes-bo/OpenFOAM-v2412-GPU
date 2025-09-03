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

#include "fileOperation.H"
#include "dummyFileOperation.H"
#include "uncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::refPtr<Foam::fileOperation> Foam::fileOperation::dummyHandlerPtr_;

Foam::refPtr<Foam::fileOperation> Foam::fileOperation::fileHandlerPtr_;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::refPtr<Foam::fileOperation> Foam::fileOperation::null()
{
    if (!dummyHandlerPtr_)
    {
        // verbose = false
        dummyHandlerPtr_.reset(new fileOperations::dummyFileOperation(false));
    }

    return dummyHandlerPtr_;
}


const Foam::fileOperation& Foam::fileOperation::fileHandler()
{
    if (!fileOperation::fileHandlerPtr_)
    {
        word handlerType(Foam::getEnv("FOAM_FILEHANDLER"));

        if (handlerType.empty())
        {
            handlerType = defaultFileHandler;
        }

        fileOperation::fileHandlerPtr_ = fileOperation::New(handlerType, true);
    }

    return *fileOperation::fileHandlerPtr_;
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(std::nullptr_t)
{
    return refPtr<fileOperation>(std::move(fileHandlerPtr_));
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(refPtr<fileOperation>& newHandler)
{
    // - do nothing if newHandler is empty. Does not delete current
    // - do nothing if newHandler is identical to current handler

    // Change ownership as atomic operations

    // If newHandler and current handler are actually identical, we
    // have a bit problem somewhere else since this means that the pointer
    // is managed is done in two places!
    // Should flag as a FatalError (in the future), but there may still be
    // some place where we would like to fake shared pointers?

    // TBD: add a flush() operation on the old handler first,
    // instead of waiting for it to be run on destruction?

    refPtr<fileOperation> old;

    if
    (
        newHandler.get() != nullptr
     && newHandler.get() != fileOperation::fileHandlerPtr_.get()
    )
    {
        old.swap(newHandler);
        old.swap(fileOperation::fileHandlerPtr_);
    }

    return old;
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(autoPtr<fileOperation>&& newHandler)
{
    // Same logic as refPtr version

    refPtr<fileOperation> old;

    if
    (
        newHandler.get() != nullptr
     && newHandler.get() != fileOperation::fileHandlerPtr_.get()
    )
    {
        old.reset(newHandler.release());
        old.swap(fileOperation::fileHandlerPtr_);
    }

    return old;
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(refPtr<fileOperation>&& newHandler)
{
    return fileOperation::fileHandler(newHandler);
}


Foam::autoPtr<Foam::fileOperation> Foam::fileOperation::NewUncollated()
{
    return autoPtr<fileOperation>
    (
        new fileOperations::uncollatedFileOperation(false)  // verbose = false
    );
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fileOperation>
Foam::fileOperation::New
(
    const word& handlerType,
    bool verbose
)
{
    if (handlerType.empty())
    {
        if (fileOperation::defaultFileHandler.empty())
        {
            FatalErrorInFunction
                << "Default file-handler name is undefined" << nl
                << abort(FatalError);
        }

        // Forward to self
        return fileOperation::New(fileOperation::defaultFileHandler, verbose);
    }

    DebugInFunction
        << "Constructing fileHandler: " << handlerType << endl;

    auto* ctorPtr = wordConstructorTable(handlerType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "fileHandler",
            handlerType,
            *wordConstructorTablePtr_
        ) << abort(FatalError);
    }

    return autoPtr<fileOperation>(ctorPtr(verbose));
}


Foam::autoPtr<Foam::fileOperation>
Foam::fileOperation::New
(
    const word& handlerType,
    const Tuple2<label, labelList>& commAndIORanks,
    const bool distributedRoots,
    bool verbose
)
{
    if (handlerType.empty())
    {
        if (fileOperation::defaultFileHandler.empty())
        {
            FatalErrorInFunction
                << "defaultFileHandler name is undefined" << nl
                << abort(FatalError);
        }

        // Forward to self
        return fileOperation::New
        (
            fileOperation::defaultFileHandler,
            commAndIORanks,
            distributedRoots,
            verbose
        );
    }

    DebugInFunction
        << "Constructing fileHandler: " << handlerType << endl;

    auto* ctorPtr = commConstructorTable(handlerType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "fileHandler",
            handlerType,
            *commConstructorTablePtr_
        ) << abort(FatalError);
    }

    return autoPtr<fileOperation>
    (
        ctorPtr(commAndIORanks, distributedRoots, verbose)
    );
}



// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// From boolUList/bitSet to list of labels,
// always include rank 0 and constrain by numProcs
template<class BoolListType>
static labelList getSelectedProcs(const BoolListType& useProc)
{
    labelList ranks;

    if
    (
        UPstream::master(UPstream::worldComm)
     || useProc.test(UPstream::myProcNo(UPstream::worldComm))
    )
    {
        DynamicList<label> subProcs(UPstream::nProcs(UPstream::worldComm));

        for (const int proci : UPstream::allProcs(UPstream::worldComm))
        {
            // Always include the master rank
            if (!proci || useProc.test(proci))
            {
                subProcs.push_back(proci);
            }
        }

        ranks.transfer(subProcs);
    }

    return ranks;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fileOperation>
Foam::fileOperation::New_impl
(
    const fileOperation& origHandler,
    const labelUList& subProcs,  // in worldComm
    bool verbose
)
{
    autoPtr<fileOperation> newHandler;

    // NB: input must include master!

    const label myProci = UPstream::myProcNo(UPstream::worldComm);
    const label numProcs = UPstream::nProcs(UPstream::worldComm);

    if (subProcs.contains(myProci))
    {
        // Retaining the original IO ranks if possible

        // Retain the original IO ranks that coincide with the new subset.
        // This may still need more attention...

        const labelUList& origIOranks = origHandler.ioRanks();
        DynamicList<label> subIORanks(origIOranks.size());

        for (const label proci : subProcs)
        {
            if (origIOranks.contains(proci))
            {
                subIORanks.push_back(proci);
            }
        }

        // Default starting point
        Tuple2<label, labelList> commAndIORanks
        (
            UPstream::worldComm,
            subIORanks
        );

        // TBD: special handling for uncollated
        // if (origHandler.comm() == UPstream::commSelf())
        // {
        //     commAndIORanks.first() = UPstream::commSelf();
        // }

        const bool hasIOranks = (commAndIORanks.second().size() > 1);

        if
        (
            UPstream::parRun()
         && (hasIOranks || (subProcs.size() != numProcs))
        )
        {
            // Without any IO range, restrict to overall proc range
            // since we don't necessarily trust the input...
            labelRange siblingRange(numProcs);

            if (hasIOranks)
            {
                // Multiple masters: ranks included in my IO range
                siblingRange = fileOperation::subRanks(commAndIORanks.second());
            }

            // Restrict to siblings within the IO range or proc range
            labelList siblings;
            if (siblingRange.size())
            {
                auto& dynSiblings = subIORanks;
                dynSiblings.clear();

                for (const label proci : subProcs)
                {
                    if (siblingRange.contains(proci))
                    {
                        dynSiblings.push_back(proci);
                    }
                }

                siblings.transfer(dynSiblings);
            }

            // Warning: MS-MPI currently uses MPI_Comm_create() instead of
            // MPI_Comm_create_group() so it will block there!

            commAndIORanks.first() = UPstream::allocateCommunicator
            (
                UPstream::worldComm,
                siblings
            );
        }


        // Allocate new handler with same type and similar IO ranks
        // but with different sub-ranks (and communicator)


        // Temporarily override world comm to get through the initialisation
        // (e.g. fileOperation::printRanks() which uses worldComm)
        const label oldWorldComm = UPstream::commWorld(commAndIORanks.first());

        newHandler = fileOperation::New
        (
            origHandler.type(),
            commAndIORanks,
            origHandler.distributed(),
            verbose
        );

        UPstream::commWorld(oldWorldComm);

        if (newHandler)
        {
            // Make sure that the output format uses the correct number of
            // 'active' ranks (instead of that of the origHandler)
            newHandler->nProcs(subProcs.size());
            newHandler->storeComm();
        }
    }

    return newHandler;
}


Foam::autoPtr<Foam::fileOperation>
Foam::fileOperation::New
(
    const fileOperation& origHandler,
    const boolUList& useProc,  // in worldComm
    bool verbose
)
{
    const labelList subProcs = getSelectedProcs(useProc);

    return New_impl(origHandler, subProcs, verbose);
}


Foam::autoPtr<Foam::fileOperation>
Foam::fileOperation::New
(
    const fileOperation& origHandler,
    const bitSet& useProc,  // in worldComm
    bool verbose
)
{
    const labelList subProcs = getSelectedProcs(useProc);

    return New_impl(origHandler, subProcs, verbose);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fileOperation>
Foam::fileHandler(autoPtr<fileOperation>&& newHandler)
{
    refPtr<fileOperation> oldHandler
    (
        fileOperation::fileHandler(std::move(newHandler))
    );

    autoPtr<fileOperation> old;

    // Can return as autoPtr if handler was also a pointer (not a reference)
    if (oldHandler.is_pointer())
    {
        old.reset(oldHandler.release());
    }

    return old;
}


// ************************************************************************* //
