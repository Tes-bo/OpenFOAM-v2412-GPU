/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

Class
    Foam::OFstreamCollator

Description
    Threaded file writer.

    Collects all data from all processors and writes as single
    'decomposedBlockData' file. The operation is determined by the
    buffer size (maxThreadFileBufferSize setting):
    - local size of data is larger than buffer: receive and write processor
    by processor (i.e. 'scheduled'). Does not use a thread, no file size
    limit.
    - total size of data is larger than buffer (but local is not):
    thread does all the collecting and writing of the processors. No file
    size limit.
    - total size of data is less than buffer:
    collecting is done locally; the thread only does the writing
    (since the data has already been collected)

SourceFiles
    OFstreamCollator.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OFstreamCollator_H
#define Foam_OFstreamCollator_H

#include <thread>
#include <mutex>
#include "IOstream.H"
#include "labelList.H"
#include "FIFOStack.H"
#include "SubList.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class OFstreamCollator Declaration
\*---------------------------------------------------------------------------*/

class OFstreamCollator
{
    // Private Class

        struct writeData
        {
            const label comm_;
            const word objectType_;
            const fileName pathName_;
            const string data_;
            const labelList sizes_;
            PtrList<List<char>> slaveData_;
            const IOstreamOption streamOpt_;
            IOstreamOption::atomicType atomic_;
            IOstreamOption::appendType append_;
            const dictionary headerEntries_;

            writeData
            (
                const label comm,
                const word& objectType,
                const fileName& pathName,
                const string& data,
                const labelList& sizes,
                IOstreamOption streamOpt,
                IOstreamOption::atomicType atomic,
                IOstreamOption::appendType append,
                const dictionary& headerEntries
            )
            :
                comm_(comm),
                objectType_(objectType),
                pathName_(pathName),
                data_(data),
                sizes_(sizes),
                slaveData_(),
                streamOpt_(streamOpt),
                atomic_(atomic),
                append_(append),
                headerEntries_(headerEntries)
            {}

            //- The (approximate) size of master + any optional slave data
            off_t size() const
            {
                off_t totalSize = data_.size();
                forAll(slaveData_, i)
                {
                    if (slaveData_.set(i))
                    {
                        totalSize += slaveData_[i].size();
                    }
                }
                return totalSize;
            }
        };


    // Private Data

        //- Total amount of storage to use for object stack below
        const off_t maxBufferSize_;

        mutable std::mutex mutex_;

        std::unique_ptr<std::thread> thread_;

        //- Stack of files to write + contents
        FIFOStack<writeData*> objects_;

        //- Whether thread is running (and not exited)
        bool threadRunning_;

        //- Communicator to use for all parallel ops (in simulation thread)
        label localComm_;

        //- Communicator to use for all parallel ops (in write thread)
        label threadComm_;


    // Private Member Functions

        //- Write actual file
        static bool writeFile
        (
            const label comm,
            const word& objectType,
            const fileName& fName,
            const string& masterData,
            const labelUList& recvSizes,
            const UPtrList<SubList<char>>& slaveData,
            IOstreamOption streamOpt,
            IOstreamOption::atomicType atomic,
            IOstreamOption::appendType append,
            const dictionary& headerEntries
        );

        //- Write all files in stack
        static void* writeAll(void *threadarg);

        //- Wait for total size of objects_ (master + optional slave data)
        //  to be wantedSize less than overall maxBufferSize.
        void waitForBufferSpace(const off_t wantedSize) const;


public:

    // Declare name of the class and its debug switch
    TypeName("OFstreamCollator");


    // Constructors

        //- Construct from buffer size. 0 = do not use thread
        explicit OFstreamCollator(const off_t maxBufferSize);

        //- Construct from buffer size (0 = do not use thread)
        //- and specified communicator
        OFstreamCollator(const off_t maxBufferSize, const label comm);


    //- Destructor
    virtual ~OFstreamCollator();


    // Member Functions

        //- Write file with contents.
        //  Blocks until writethread has space available
        //  (total file sizes < maxBufferSize)
        bool write
        (
            const word& objectType,
            const fileName&,
            const string& data,
            IOstreamOption streamOpt,
            IOstreamOption::atomicType atomic,
            IOstreamOption::appendType append,
            const bool useThread = true,
            const dictionary& headerEntries = dictionary::null
        );

        //- Wait for all thread actions to have finished
        void waitAll();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
