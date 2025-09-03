/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "collatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "threadedCollatedOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "masterOFstream.H"
#include "OFstream.H"
#include "foamVersion.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(collatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        collatedFileOperation,
        word
    );
    addToRunTimeSelectionTable
    (
        fileOperation,
        collatedFileOperation,
        comm
    );

    float collatedFileOperation::maxThreadFileBufferSize
    (
        debug::floatOptimisationSwitch("maxThreadFileBufferSize", 0)
    );
    registerOptSwitch
    (
        "maxThreadFileBufferSize",
        float,
        collatedFileOperation::maxThreadFileBufferSize
    );

    // Threaded MPI: depending on buffering
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        fileOperationInitialise_collated,
        word,
        collated
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileOperations::collatedFileOperation::printBanner
(
    const bool withRanks
) const
{
    DetailInfo
        << "I/O    : " << this->type();

    if (mag(maxThreadFileBufferSize) > 1)
    {
        // FUTURE: deprecate or remove threading?
        DetailInfo
            << " [threaded] (maxThreadFileBufferSize = "
            << maxThreadFileBufferSize << ")." << nl
            << "         Requires buffer large enough to collect all data"
               " or MPI thread support." << nl
            << "         To avoid MPI threading [slow], set"
               " (maxThreadFileBufferSize = 0) in" << nl
            << "         OpenFOAM etc/controlDict" << endl;
    }
    else
    {
        DetailInfo
            << " [unthreaded] (maxThreadFileBufferSize = 0)." << nl;

        if (mag(maxMasterFileBufferSize) < 1)
        {
            DetailInfo
                << "         With scheduled transfer" << nl;
        }
        else if (maxMasterFileBufferSize >= 1)
        {
            DetailInfo
                << "         With non-blocking transfer,"
                   " buffer-size = " << maxMasterFileBufferSize << nl;
        }
        else
        {
            DetailInfo
                << "         With non-blocking transfer,"
                   " minimal buffer size" << nl;
        }
    }

    if (withRanks)
    {
        fileOperation::printRanks();
    }

    //- fileModificationChecking already set by base class (masterUncollated)
    // if (IOobject::fileModificationChecking == IOobject::timeStampMaster)
    // {
    //     WarningInFunction
    //         << "Resetting fileModificationChecking to timeStamp" << endl;
    // }
    // else if (IOobject::fileModificationChecking == IOobject::inotifyMaster)
    // {
    //     WarningInFunction
    //         << "Resetting fileModificationChecking to inotify" << endl;
    // }
}


bool Foam::fileOperations::collatedFileOperation::appendObject
(
    const regIOobject& io,
    const fileName& pathName,
    IOstreamOption streamOpt
) const
{
    // Append to processorsNN/ file

    const label proci = detectProcessorPath(io.objectPath());

    if (debug)
    {
        Pout<< "collatedFileOperation::writeObject :"
            << " For local object : " << io.name()
            << " appending processor " << proci
            << " data to " << pathName << endl;
    }
    if (proci == -1)
    {
        FatalErrorInFunction
            << "Invalid processor path: " << pathName
            << exit(FatalError);
    }

    const bool isIOmaster = fileOperation::isIOrank(proci);

    // Update meta-data for current state
    if (isIOmaster)
    {
        const_cast<regIOobject&>(io).updateMetaData();
    }

    // Note: cannot do append + compression. This is a limitation
    // of ogzstream (or rather most compressed formats)
    //
    // File should always be created as non-atomic
    // (consistency between append/non-append)

    OFstream os
    (
        pathName,
        // UNCOMPRESSED (binary only)
        IOstreamOption(IOstreamOption::BINARY, streamOpt.version()),
        // Append on sub-ranks
        (isIOmaster ? IOstreamOption::NO_APPEND : IOstreamOption::APPEND_APP)
    );

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Cannot open for appending"
            << exit(FatalIOError);
    }

    if (isIOmaster)
    {
        decomposedBlockData::writeHeader(os, streamOpt, io);
    }

    std::streamoff blockOffset = decomposedBlockData::writeBlockEntry
    (
        os,
        streamOpt,
        io,
        proci,
        isIOmaster  // With FoamFile header on master
    );

    return (blockOffset >= 0) && os.good();
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

void Foam::fileOperations::collatedFileOperation::init(bool verbose)
{
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (verbose)
    {
        this->printBanner(ioRanks_.size());
    }
}


Foam::fileOperations::collatedFileOperation::collatedFileOperation
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
    managedComm_(getManagedComm(comm_)),  // Possibly locally allocated
    writer_(mag(maxThreadFileBufferSize), comm_)
{
    init(verbose);
}


Foam::fileOperations::collatedFileOperation::collatedFileOperation
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
    managedComm_(-1),  // Externally managed
    writer_(mag(maxThreadFileBufferSize), comm_)
{
    init(verbose);
}


void Foam::fileOperations::collatedFileOperation::storeComm() const
{
    // From externally -> locally managed
    managedComm_ = getManagedComm(comm_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::~collatedFileOperation()
{
    // Wait for any outstanding file operations
    flush();

    UPstream::freeCommunicator(managedComm_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::collatedFileOperation::objectPath
(
    const IOobject& io,
    const word& typeName
) const
{
    // Replacement for objectPath
    if (io.time().processorCase())
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::PROCOBJECT,
            "dummy",        // not used for processorsobject
            io.instance()
        );
    }
    else
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::OBJECT,
            word::null,
            io.instance()
        );
    }
}


bool Foam::fileOperations::collatedFileOperation::writeObject
(
    const regIOobject& io,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    const Time& tm = io.time();
    const fileName& inst = io.instance();

    // Update meta-data for current state
    const_cast<regIOobject&>(io).updateMetaData();

    if (inst.isAbsolute() || !tm.processorCase())
    {
        // Note: delay mkdir to masterOFstream so it does not get created
        //       if not needed (e.g. when running distributed)

        const fileName pathName(io.objectPath());

        if (debug)
        {
            Pout<< "collatedFileOperation::writeObject :"
                << " For object : " << io.name()
                << " falling back to master-only output to " << io.path()
                << endl;
        }

        // Note: currently still NON_ATOMIC (Dec-2022)
        masterOFstream os
        (
            comm_,
            pathName,
            streamOpt,
            IOstreamOption::NO_APPEND,
            writeOnProc
        );

        // If any of these fail, return
        // (leave error handling to Ostream class)

        const bool ok =
        (
            os.good()
         && io.writeHeader(os)
         && io.writeData(os)
        );

        if (ok)
        {
            IOobject::writeEndDivider(os);
        }

        return ok;
    }
    else
    {
        // Construct the equivalent processors/ directory
        const fileName path(processorsPath(io, inst, processorsDir(io)));

        // Note: delay mkdir to masterOFstream so it does not get created
        //       if not needed (e.g. when running distributed)

        const fileName pathName(path/io.name());

        if (io.global() || io.globalObject())
        {
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For global object : " << io.name()
                    << " falling back to master-only output to " << pathName
                    << endl;
            }

            // Note: currently still NON_ATOMIC (Dec-2022)
            masterOFstream os
            (
                comm_,
                pathName,
                streamOpt,
                IOstreamOption::NO_APPEND,
                writeOnProc
            );

            // If any of these fail, return
            // (leave error handling to Ostream class)

            const bool ok =
            (
                os.good()
             && io.writeHeader(os)
             && io.writeData(os)
            );

            if (ok)
            {
                IOobject::writeEndDivider(os);
            }

            return ok;
        }
        else if (!UPstream::parRun())
        {
            // Special path for e.g. decomposePar. Append to
            // processorsDDD/ file
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " appending to " << pathName << endl;
            }

            mkDir(path);
            return appendObject(io, pathName, streamOpt);
        }
        else
        {
            // Re-check static maxThreadFileBufferSize variable to see
            // if needs to use threading
            const bool useThread = (mag(maxThreadFileBufferSize) > 1);

            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " starting collating output to " << pathName
                    << " useThread:" << useThread << endl;
            }

            if (!useThread)
            {
                writer_.waitAll();
            }

            // Note: currently still NON_ATOMIC (Dec-2022)
            threadedCollatedOFstream os
            (
                writer_,
                pathName,
                streamOpt,
                useThread
            );

            bool ok = os.good();

            if (UPstream::master(comm_))
            {
                // Suppress comment banner
                const bool old = IOobject::bannerEnabled(false);

                ok = ok && io.writeHeader(os);

                IOobject::bannerEnabled(old);

                // Additional header content
                dictionary dict;
                decomposedBlockData::writeExtraHeaderContent
                (
                    dict,
                    streamOpt,
                    io
                );
                os.setHeaderEntries(dict);
            }

            ok = ok && io.writeData(os);
            // No end divider for collated output

            return ok;
        }
    }
}

void Foam::fileOperations::collatedFileOperation::flush() const
{
    if (debug)
    {
        Pout<< "collatedFileOperation::flush : clearing and waiting for thread"
            << endl;
    }
    masterUncollatedFileOperation::flush();
    // Wait for thread to finish (note: also removes thread)
    writer_.waitAll();
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const fileName& fName
) const
{
    if (UPstream::parRun())
    {
        const List<int>& procs(UPstream::procID(comm_));

        word procDir(processorsBaseDir+Foam::name(nProcs_));

        if (procs.size() != nProcs_)
        {
            procDir +=
              + "_"
              + Foam::name(procs.first())
              + "-"
              + Foam::name(procs.last());
        }
        return procDir;
    }
    else
    {
        word procDir(processorsBaseDir+Foam::name(nProcs_));

        if (ioRanks_.size())
        {
            // Detect current processor number
            label proci = detectProcessorPath(fName);

            if (proci != -1)
            {
                // Find lowest io rank
                label minProc = 0;
                label maxProc = nProcs_-1;
                for (const label ranki : ioRanks_)
                {
                    if (ranki >= nProcs_)
                    {
                        break;
                    }
                    else if (ranki <= proci)
                    {
                        minProc = ranki;
                    }
                    else
                    {
                        maxProc = ranki-1;
                        break;
                    }
                }

                // Add range if not all processors
                if (maxProc-minProc+1 != nProcs_)
                {
                    procDir +=
                      + "_"
                      + Foam::name(minProc)
                      + "-"
                      + Foam::name(maxProc);
                }
            }
        }

        return procDir;
    }
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const IOobject& io
) const
{
    return processorsDir(io.objectPath());
}


// ************************************************************************* //
