/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "OFstreamCollator.H"
#include "OFstream.H"
#include "decomposedBlockData.H"
#include "dictionary.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstreamCollator, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::OFstreamCollator::writeFile
(
    const label comm,
    const word& objectType,
    const fileName& fName,
    const string& masterData,
    const labelUList& recvSizes,
    const UPtrList<SubList<char>>& slaveData,    // optional slave data
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const dictionary& headerEntries
)
{
    if (debug)
    {
        Pout<< "OFstreamCollator : Writing master " << label(masterData.size())
            << " bytes to " << fName << " using comm " << comm
            << " and " << slaveData.size() << " sub-ranks" << endl;

        forAll(slaveData, proci)
        {
            if (slaveData.set(proci))
            {
                Pout<< "    " << proci
                    << " size:" << slaveData[proci].size()
                    << endl;
            }
        }
    }

    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm))
    {
        Foam::mkDir(fName.path());
        osPtr.reset(new OFstream(atomic, fName, streamOpt, append));
        auto& os = *osPtr;

        if (append == IOstreamOption::NO_APPEND)
        {
            // No IOobject so cannot use IOobject::writeHeader

            // FoamFile
            decomposedBlockData::writeHeader
            (
                os,
                streamOpt,      // streamOpt for container
                objectType,
                "",             // note
                "",             // location (leave empty instead inaccurate)
                fName.name(),   // object name
                headerEntries
            );
        }
    }

    // Assuming threaded writing hides any slowness so we
    // can use scheduled communication to send the data to
    // the master processor in order. However can be unstable
    // for some mpi so default is non-blocking.
    const UPstream::commsTypes myCommsType
    (
        (
            fileOperations::masterUncollatedFileOperation::
                maxMasterFileBufferSize == 0
        )
      ? UPstream::commsTypes::scheduled
      : UPstream::commsTypes::nonBlocking
    );


    UList<char> slice
    (
        const_cast<char*>(masterData.data()),
        label(masterData.size())
    );

    List<std::streamoff> blockOffset;
    decomposedBlockData::writeBlocks
    (
        comm,
        osPtr,
        blockOffset,
        slice,
        recvSizes,
        slaveData,
        myCommsType,
        false       // do not reduce return state
    );

    if (osPtr && !osPtr->good())
    {
        FatalIOErrorInFunction(*osPtr)
            << "Failed writing to " << fName << exit(FatalIOError);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Finished writing " << masterData.size()
            << " bytes";
        if (UPstream::master(comm))
        {
            off_t sum = 0;
            for (const label recv : recvSizes)
            {
                sum += recv;
            }
            // Use std::to_string to display long int
            Pout<< " (overall " << std::to_string(sum) << ')';
        }
        Pout<< " to " << fName
            << " using comm " << comm << endl;
    }

    return true;
}


void* Foam::OFstreamCollator::writeAll(void *threadarg)
{
    OFstreamCollator& handler = *static_cast<OFstreamCollator*>(threadarg);

    // Consume stack
    while (true)
    {
        writeData* ptr = nullptr;

        {
            std::lock_guard<std::mutex> guard(handler.mutex_);
            if (handler.objects_.size())
            {
                ptr = handler.objects_.pop();
            }
        }

        if (!ptr)
        {
            break;
        }
        else
        {
            // Convert storage to pointers
            PtrList<SubList<char>> slaveData;
            if (ptr->slaveData_.size())
            {
                slaveData.resize(ptr->slaveData_.size());
                forAll(slaveData, proci)
                {
                    if (ptr->slaveData_.set(proci))
                    {
                        slaveData.set
                        (
                            proci,
                            new SubList<char>
                            (
                                ptr->slaveData_[proci],
                                ptr->sizes_[proci]
                            )
                        );
                    }
                }
            }

            bool ok = writeFile
            (
                ptr->comm_,
                ptr->objectType_,
                ptr->pathName_,
                ptr->data_,
                ptr->sizes_,
                slaveData,
                ptr->streamOpt_,
                ptr->atomic_,
                ptr->append_,
                ptr->headerEntries_
            );
            if (!ok)
            {
                FatalIOErrorInFunction(ptr->pathName_)
                    << "Failed writing " << ptr->pathName_
                    << exit(FatalIOError);
            }

            delete ptr;
        }
        //sleep(1);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Exiting write thread " << endl;
    }

    {
        std::lock_guard<std::mutex> guard(handler.mutex_);
        handler.threadRunning_ = false;
    }

    return nullptr;
}


void Foam::OFstreamCollator::waitForBufferSpace(const off_t wantedSize) const
{
    while (true)
    {
        // Count files to be written
        off_t totalSize = 0;

        {
            std::lock_guard<std::mutex> guard(mutex_);
            forAllConstIters(objects_, iter)
            {
                totalSize += iter()->size();
            }
        }

        if
        (
            totalSize == 0
         || (wantedSize >= 0 && (totalSize+wantedSize) <= maxBufferSize_)
        )
        {
            break;
        }

        if (debug)
        {
            std::lock_guard<std::mutex> guard(mutex_);
            Pout<< "OFstreamCollator : Waiting for buffer space."
                << " Currently in use:" << totalSize
                << " limit:" << maxBufferSize_
                << " files:" << objects_.size()
                << endl;
        }

        sleep(5);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstreamCollator::OFstreamCollator(const off_t maxBufferSize)
:
    maxBufferSize_(maxBufferSize),
    threadRunning_(false),
    localComm_(UPstream::worldComm),
    threadComm_
    (
        UPstream::allocateCommunicator
        (
            localComm_,
            labelRange(UPstream::nProcs(localComm_))
        )
    )
{}


Foam::OFstreamCollator::OFstreamCollator
(
    const off_t maxBufferSize,
    const label comm
)
:
    maxBufferSize_(maxBufferSize),
    threadRunning_(false),
    localComm_(comm),
    threadComm_
    (
        UPstream::allocateCommunicator
        (
            localComm_,
            labelRange(UPstream::nProcs(localComm_))
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstreamCollator::~OFstreamCollator()
{
    if (thread_)
    {
        if (debug)
        {
            Pout<< "~OFstreamCollator : Waiting for write thread" << endl;
        }
        thread_->join();
        thread_.reset(nullptr);
    }

    UPstream::freeCommunicator(threadComm_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OFstreamCollator::write
(
    const word& objectType,
    const fileName& fName,
    const string& data,
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const bool useThread,
    const dictionary& headerEntries
)
{
    // Determine (on master) sizes to receive. Note: do NOT use thread
    // communicator
    labelList recvSizes;
    decomposedBlockData::gather(localComm_, label(data.size()), recvSizes);

    off_t totalSize = 0;
    label maxLocalSize = 0;
    {
        if (UPstream::master(localComm_))
        {
            for (const label recvSize : recvSizes)
            {
                totalSize += recvSize;
                maxLocalSize = max(maxLocalSize, recvSize);
            }
        }
        Pstream::broadcasts(localComm_, totalSize, maxLocalSize);
    }

    if (!useThread || maxBufferSize_ == 0 || maxLocalSize > maxBufferSize_)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather and write of " << fName
                << " using local comm " << localComm_ << endl;
        }
        // Direct collating and writing (so master blocks until all written!)
        const PtrList<SubList<char>> dummySlaveData;
        return writeFile
        (
            localComm_,
            objectType,
            fName,
            data,
            recvSizes,
            dummySlaveData,
            streamOpt,
            atomic,
            append,
            headerEntries
        );
    }
    else if (totalSize <= maxBufferSize_)
    {
        // Total size can be stored locally so receive all data now and only
        // do the writing in the thread

        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather; thread write of "
                << fName << endl;
        }

        if (Pstream::master(localComm_))
        {
            waitForBufferSpace(totalSize);
        }


        // Receive in chunks of labelMax (2^31-1) since this is the maximum
        // size that a List can be

        autoPtr<writeData> fileAndDataPtr
        (
            new writeData
            (
                threadComm_,        // Note: comm not actually used anymore
                objectType,
                fName,
                (
                    Pstream::master(localComm_)
                  ? data            // Only used on master
                  : string::null
                ),
                recvSizes,
                streamOpt,
                atomic,
                append,
                headerEntries
            )
        );
        writeData& fileAndData = fileAndDataPtr();

        PtrList<List<char>>& slaveData = fileAndData.slaveData_;

        UList<char> slice(const_cast<char*>(data.data()), label(data.size()));

        slaveData.setSize(recvSizes.size());

        // Gather all data onto master. Is done in local communicator since
        // not in write thread. Note that we do not store in contiguous
        // buffer since that would limit to 2G chars.
        const label startOfRequests = UPstream::nRequests();
        if (Pstream::master(localComm_))
        {
            for (label proci = 1; proci < slaveData.size(); proci++)
            {
                slaveData.set(proci, new List<char>(recvSizes[proci]));
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    slaveData[proci].data(),
                    slaveData[proci].size_bytes(),
                    Pstream::msgType(),
                    localComm_
                );
            }
        }
        else
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    0,
                    slice.cdata(),
                    slice.size_bytes(),
                    Pstream::msgType(),
                    localComm_
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << 0 << " nBytes:"
                    << label(slice.size_bytes())
                    << Foam::abort(FatalError);
            }
        }
        UPstream::waitRequests(startOfRequests);

        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Append to thread buffer
            objects_.push(fileAndDataPtr.ptr());

            // Start thread if not running
            if (!threadRunning_)
            {
                if (thread_)
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_->join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread"
                        << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }
    else
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : thread gather and write of " << fName
                << " using communicator " << threadComm_ << endl;
        }

        if (!UPstream::haveThreads())
        {
            FatalErrorInFunction
                << "mpi does not seem to have thread support."
                << " Make sure to set buffer size 'maxThreadFileBufferSize'"
                << " to at least " << totalSize
                << " to be able to do the collating before threading."
                << exit(FatalError);
        }

        if (Pstream::master(localComm_))
        {
            waitForBufferSpace(data.size());
        }

        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Push all file info on buffer. Note that no slave data provided
            // so it will trigger communication inside the thread
            objects_.push
            (
                new writeData
                (
                    threadComm_,
                    objectType,
                    fName,
                    data,
                    recvSizes,
                    streamOpt,
                    atomic,
                    append,
                    headerEntries
                )
            );

            if (!threadRunning_)
            {
                if (thread_)
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_->join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread" << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }
}


void Foam::OFstreamCollator::waitAll()
{
    // Wait for all buffer space to be available i.e. wait for all jobs
    // to finish
    if (Pstream::master(localComm_))
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : waiting for thread to have consumed all"
                << endl;
        }
        waitForBufferSpace(-1);
    }
}


// ************************************************************************* //
