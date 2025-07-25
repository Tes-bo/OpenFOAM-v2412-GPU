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

#include "decomposedBlockData.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamBuffers.H"
#include "Fstream.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "SubList.H"
#include "charList.H"
#include "labelPair.H"
#include "masterUncollatedFileOperation.H"
#include "SpanStream.H"
#include "StringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::decomposedBlockData::isCollatedType
(
    const word& objectType
)
{
    return
    (
        objectType == decomposedBlockData::typeName
    );
}


bool Foam::decomposedBlockData::isCollatedType
(
    const IOobject& io
)
{
    return decomposedBlockData::isCollatedType(io.headerClassName());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decomposedBlockData::decomposedBlockData
(
    const label comm,
    const IOobject& io,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(comm),
    contentData_()
{
    // Temporary warning
    if (readOpt() == IOobjectOption::READ_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with READ_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        read();
    }
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    const label blocki,
    const char* str,
    const size_t len
)
{
    // Offset to the beginning of this output
    // This should generally be OK for non-compressed streams
    // (eg, std::ofstream)

    std::streamoff blockOffset = os.stdStream().tellp();

    const word procName("processor" + Foam::name(blocki));

    // Write as commented content
    // ----------------
    // // processorN
    // NCHARS
    // (...)
    // ----------------
    {
        os << nl << "// " << procName << nl;

        if (str && len > 0)
        {
            // Special treatment for char data (binary I/O only)
            const auto oldFmt = os.format(IOstreamOption::BINARY);

            os << label(len) << nl;
            os.write(str, len);
            os << nl;

            os.format(oldFmt);
        }
        else
        {
            os << label(0) << nl;
        }
    }

    // Write as primitiveEntry
    // {
    //     // Like writeKeyword()
    //     os << nl << procName << nl;
    //
    //     if (str && len > 0)
    //     {
    //         // Special treatment for char data (binary I/O only)
    //         const auto oldFmt = os.format(IOstreamOption::BINARY);
    //
    //         os << label(len) << nl;
    //         os.write(str, len);
    //         os << nl;
    //
    //         os.format(oldFmt);
    //     }
    //     else
    //     {
    //         os << label(0) << nl;
    //     }
    //
    //     os.endEntry();
    // }

    return blockOffset;
}


bool Foam::decomposedBlockData::readBlockEntry
(
    Istream& is,
    List<char>& charData
)
{
    // Handle any of these:

    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;
    // 3.  processorN  NCHARS (...) ;

    is.fatalCheck(FUNCTION_NAME);
    token tok(is);
    is.fatalCheck(FUNCTION_NAME);

    // Dictionary format has primitiveEntry keyword:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (!isDictFormat && tok.good())
    {
        is.putBack(tok);
    }
    charData.readList(is);

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return true;
}


bool Foam::decomposedBlockData::skipBlockEntry(Istream& is)
{
    // As per readBlockEntry but seeks instead of reading.
    // Internals like charList::readList - ie, always binary

    // Handle any of these:
    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;
    // 3.  processorN  NCHARS (...) ;

    if (!is.good()) return false;
    token tok(is);
    if (!is.good()) return false;

    // Dictionary format has primitiveEntry keyword:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (isDictFormat)
    {
        is >> tok;
        if (!is.good()) return false;
    }


    bool handled = false;

    // Like charList::readList
    if (tok.isCompound())
    {
        handled = true;
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..) or just a plain '0'

        const label len = tok.labelToken();

        // Special treatment for char data (binary I/O only)
        const auto oldFmt = is.format(IOstreamOption::BINARY);

        if (len)
        {
            // read(...) includes surrounding start/end delimiters.

            // Note: nullptr to ignore instead of reading
            is.read(nullptr, std::streamsize(len));
        }
        is.format(oldFmt);

        handled = true;
    }
    else
    {
        // Incorrect token
        return false;
    }

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return handled;
}


Foam::label Foam::decomposedBlockData::getNumBlocks
(
    Istream& is,
    const label maxNumBlocks
)
{
    label nBlocks = 0;

    // Handle OpenFOAM header if it is the first entry
    if (is.good())
    {
        token tok(is);

        if (is.good() && tok.isWord("FoamFile"))
        {
            dictionary headerDict(is);  // Read sub-dictionary content

            if (headerDict.readIfPresent("version", tok))
            {
                is.version(tok);
            }

            word formatName;
            if (headerDict.readIfPresent("format", formatName))
            {
                is.format(formatName);
            }

            //// Obtain number of blocks directly
            ///  This may not be reliable...
            //if (headerDict.readIfPresent("blocks", nBlocks))
            //{
            //    return nBlocks;
            //}
        }
        else if (tok.good())
        {
            is.putBack(tok);
        }
    }

    while (is.good() && skipBlockEntry(is))
    {
        ++nBlocks;

        if (maxNumBlocks == nBlocks)
        {
            break;
        }
    }

    return nBlocks;
}


bool Foam::decomposedBlockData::hasBlock(Istream& is, const label blockNumber)
{
    return
    (
        blockNumber >= 0
     && (blockNumber < getNumBlocks(is, blockNumber+1))
    );
}


std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    IOstreamOption streamOptData,
    const regIOobject& io,
    const label blocki,
    const bool withLocalHeader
)
{
    // String of all data to write
    string contentChars;
    {
        OStringStream buf(streamOptData);

        bool ok = true;

        // Generate FoamFile header on master, without comment banner
        if (withLocalHeader)
        {
            const bool old = IOobject::bannerEnabled(false);

            ok = io.writeHeader(buf);

            IOobject::bannerEnabled(old);
        }

        // Write the data to the Ostream
        ok = ok && io.writeData(buf);

        if (!ok)
        {
            return std::streamoff(-1);
        }

        contentChars = buf.str();
    }

    return decomposedBlockData::writeBlockEntry(os, blocki, contentChars);
}


Foam::autoPtr<Foam::ISstream>
Foam::decomposedBlockData::readBlock
(
    const label blocki,
    ISstream& is,
    IOobject& headerIO
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlock:"
            << " stream:" << is.name() << " attempt to read block " << blocki
            << endl;
    }

    // Extracted header information
    IOstreamOption streamOptData;
    unsigned labelWidth = is.labelByteSize();
    unsigned scalarWidth = is.scalarByteSize();

    autoPtr<ISstream> realIsPtr;

    // Read master for header
    List<char> data;
    decomposedBlockData::readBlockEntry(is, data);

    if (blocki == 0)
    {
        realIsPtr.reset(new ICharStream(std::move(data)));
        realIsPtr->name() = is.name();

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*realIsPtr))
            {
                FatalIOErrorInFunction(*realIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }
    }
    else
    {
        {
            // Read header from first block
            ISpanStream headerStream(data);
            if (!headerIO.readHeader(headerStream))
            {
                FatalIOErrorInFunction(headerStream)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
            streamOptData = static_cast<IOstreamOption>(headerStream);
            labelWidth = headerStream.labelByteSize();
            scalarWidth = headerStream.scalarByteSize();
        }

        for (label i = 1; i < blocki+1; i++)
        {
            // Read and discard data, only retain the last one
            decomposedBlockData::readBlockEntry(is, data);
        }
        realIsPtr.reset(new ICharStream(std::move(data)));
        realIsPtr->name() = is.name();

        // Apply stream settings
        realIsPtr().format(streamOptData.format());
        realIsPtr().version(streamOptData.version());
        realIsPtr().setLabelByteSize(labelWidth);
        realIsPtr().setScalarByteSize(scalarWidth);
    }

    return realIsPtr;
}


bool Foam::decomposedBlockData::readBlocks
(
    const label comm,
    autoPtr<ISstream>& isPtr,
    List<char>& data,
    const UPstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "invalid")
            << " commsType:" << Pstream::commsTypeNames[commsType]
            << " comm:" << comm << endl;
    }

    bool ok = false;

    if (UPstream::master(comm))
    {
        auto& is = *isPtr;
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, data);
    }

    if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                List<char> elems;
                decomposedBlockData::readBlockEntry(is, elems);

                OPstream os
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    0,
                    UPstream::msgType(),
                    comm
                );
                os << elems;
            }

            ok = is.good();
        }
        else
        {
            IPstream is
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                0,
                UPstream::msgType(),
                comm
            );
            is >> data;
        }
    }
    else
    {
        PstreamBuffers pBufs(comm);

        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                List<char> elems;
                decomposedBlockData::readBlockEntry(is, elems);

                UOPstream os(proci, pBufs);
                os << elems;
            }
        }

        pBufs.finishedScatters();

        if (!UPstream::master(comm))
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> data;
        }
    }

    Pstream::broadcast(ok, comm);

    return ok;
}


Foam::autoPtr<Foam::ISstream> Foam::decomposedBlockData::readBlocks
(
    const label comm,
    const fileName& fName,
    autoPtr<ISstream>& isPtr,
    IOobject& headerIO,
    const UPstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "invalid")
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    bool ok = false;
    List<char> data;
    autoPtr<ISstream> realIsPtr;

    if (UPstream::master(comm))
    {
        auto& is = *isPtr;
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, data);

        realIsPtr.reset(new ICharStream(std::move(data)));
        realIsPtr->name() = fName;

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*realIsPtr))
            {
                FatalIOErrorInFunction(*realIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }
    }

    if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                decomposedBlockData::readBlockEntry(is, data);

                OPstream os
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    0,
                    UPstream::msgType(),
                    comm
                );
                os << data;
            }

            ok = is.good();
        }
        else
        {
            IPstream is
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                0,
                UPstream::msgType(),
                comm
            );
            is >> data;

            realIsPtr.reset(new ICharStream(std::move(data)));
            realIsPtr->name() = fName;
        }
    }
    else
    {
        PstreamBuffers pBufs(comm);

        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                List<char> elems;
                decomposedBlockData::readBlockEntry(is, elems);

                UOPstream os(proci, pBufs);
                os << elems;
            }

            ok = is.good();
        }

        pBufs.finishedScatters();

        if (!UPstream::master(comm))
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> data;

            realIsPtr.reset(new ICharStream(std::move(data)));
            realIsPtr->name() = fName;
        }
    }

    Pstream::broadcast(ok, comm);

    // Broadcast master header info,
    // set stream properties from realIsPtr on master

    int verValue;
    int fmtValue;
    unsigned labelWidth;
    unsigned scalarWidth;
    word headerName(headerIO.name());

    if (UPstream::master(comm))
    {
        verValue = realIsPtr().version().canonical();
        fmtValue = static_cast<int>(realIsPtr().format());
        labelWidth = realIsPtr().labelByteSize();
        scalarWidth = realIsPtr().scalarByteSize();
    }

    Pstream::broadcasts
    (
        UPstream::worldComm,   // Future? comm,
        verValue,
        fmtValue,
        labelWidth,
        scalarWidth,
        headerName,
        headerIO.headerClassName(),
        headerIO.note()
        // Unneeded: headerIO.instance()
        // Unneeded: headerIO.local()
    );

    realIsPtr().version(IOstreamOption::versionNumber::canonical(verValue));
    realIsPtr().format(IOstreamOption::streamFormat(fmtValue));
    realIsPtr().setLabelByteSize(labelWidth);
    realIsPtr().setScalarByteSize(scalarWidth);

    headerIO.rename(headerName);

    return realIsPtr;
}


void Foam::decomposedBlockData::gather
(
    const label comm,
    const label data,
    labelList& datas
)
{
    const label nProcs = UPstream::nProcs(comm);
    datas.resize(nProcs);

    char* data0Ptr = datas.data_bytes();

    List<int> recvOffsets;
    List<int> recvSizes;
    if (UPstream::master(comm))
    {
        recvOffsets.setSize(nProcs);
        forAll(recvOffsets, proci)
        {
            // Note: truncating long int to int since UPstream::gather limited
            // to ints
            recvOffsets[proci] =
                int(reinterpret_cast<char*>(&datas[proci]) - data0Ptr);
        }
        recvSizes.setSize(nProcs, sizeof(label));
    }

    UPstream::gather
    (
        reinterpret_cast<const char*>(&data),
        sizeof(label),
        data0Ptr,
        recvSizes,
        recvOffsets,
        comm
    );
}


void Foam::decomposedBlockData::gatherSlaveData
(
    const label comm,
    const UList<char>& data,
    const labelUList& recvSizes,

    const labelRange& fromProcs,

    List<int>& sliceOffsets,
    DynamicList<char>& recvData
)
{
    const label myProci = UPstream::myProcNo(comm);
    const label numProcs = UPstream::nProcs(comm);

    int nSendBytes = 0;
    recvData.clear();

    // Calculate master data
    List<int> sliceSizes;
    if (UPstream::master(comm))
    {
        sliceSizes.resize(numProcs, 0);
        sliceOffsets.resize(numProcs+1, 0);

        // Offset 1 beyond the end of the range
        const label endProci = fromProcs.end_value();

        int totalSize = 0;
        for (const label proci : fromProcs)
        {
            sliceSizes[proci] = int(recvSizes[proci]);
            sliceOffsets[proci] = totalSize;
            totalSize += sliceSizes[proci];
        }

        sliceOffsets[endProci] = totalSize;
        recvData.resize(totalSize);
    }
    else if (fromProcs.contains(myProci))
    {
        // Note: UPstream::gather limited to int
        nSendBytes = int(data.size_bytes());
    }

    UPstream::gather
    (
        data.cdata(),
        nSendBytes,

        recvData.data(),
        sliceSizes,
        sliceOffsets,
        comm
    );
}


Foam::label Foam::decomposedBlockData::calcNumProcs
(
    const label comm,
    const off_t maxBufferSize,
    const labelUList& recvSizes,
    const label startProci
)
{
    const label nProcs = UPstream::nProcs(comm);

    label nSendProcs = 0;
    if (UPstream::master(comm))
    {
        off_t totalSize = recvSizes[startProci];
        label proci = startProci+1;
        while (proci < nProcs && (totalSize+recvSizes[proci] < maxBufferSize))
        {
            totalSize += recvSizes[proci];
            proci++;
        }

        nSendProcs = proci-startProci;
    }

    Pstream::broadcast(nSendProcs, comm);

    return nSendProcs;
}


bool Foam::decomposedBlockData::writeBlocks
(
    const label comm,
    autoPtr<OSstream>& osPtr,
    List<std::streamoff>& blockOffset,
    const UList<char>& masterData,

    const labelUList& recvSizes,
    const UPtrList<SubList<char>>& slaveData,

    const UPstream::commsTypes commsType,
    const bool syncReturnState
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::writeBlocks:"
            << " stream:" << (osPtr ? osPtr->name() : "none")
            << " data:" << masterData.size()
            << " (master only) slaveData:" << slaveData.size()
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    const label nProcs = UPstream::nProcs(comm);

    bool ok = true;

    // Write master data
    if (UPstream::master(comm))
    {
        blockOffset.resize(nProcs);

        OSstream& os = osPtr();

        blockOffset[UPstream::masterNo()] =
            decomposedBlockData::writeBlockEntry
            (
                os,
                UPstream::masterNo(),
                masterData
            );

        ok = os.good();
    }

    if (slaveData.size())
    {
        // Already have gathered the slave data.

        if (UPstream::master(comm))
        {
            // Master data already written ...
            OSstream& os = osPtr();

            // Write slaves
            for (label proci = 1; proci < nProcs; ++proci)
            {
                blockOffset[proci] =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        slaveData[proci]
                    );
            }

            ok = os.good();
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            // Master data already written ...
            OSstream& os = osPtr();

            // Receive and write slaves
            label maxNonLocalSize = 0;
            for (label proci = 1; proci < nProcs; ++proci)
            {
                maxNonLocalSize = max(maxNonLocalSize, recvSizes[proci]);
            }

            DynamicList<char> elems(maxNonLocalSize);
            for (label proci = 1; proci < nProcs; ++proci)
            {
                elems.resize_nocopy(recvSizes[proci]);
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    elems.data(),
                    elems.size_bytes(),
                    UPstream::msgType(),
                    comm
                );

                blockOffset[proci] =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        elems
                    );
            }

            ok = os.good();
        }
        else
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                masterData.cdata(),
                masterData.size_bytes(),
                UPstream::msgType(),
                comm
            );
        }
    }
    else
    {
        // Master data already written ...

        // Find out how many processor can be received into
        // maxMasterFileBufferSize

        // Starting slave processor and number of processors
        label startProc = 1;
        label nSendProcs = nProcs-1;

        DynamicList<char> recvData;

        while (nSendProcs > 0 && startProc < nProcs)
        {
            nSendProcs = calcNumProcs
            (
                comm,
                off_t
                (
                    fileOperations::masterUncollatedFileOperation::
                    maxMasterFileBufferSize
                ),
                recvSizes,
                startProc
            );

            if (nSendProcs == 0)
            {
                break;
            }


            // Gather data from (a slice of) the slaves
            labelRange fromProcs(startProc, nSendProcs);

            List<int> sliceOffsets;
            gatherSlaveData
            (
                comm,
                masterData,
                recvSizes,

                fromProcs,

                sliceOffsets,
                recvData
            );

            if (UPstream::master(comm))
            {
                OSstream& os = osPtr();

                // Write received data
                for (const label proci : fromProcs)
                {
                    SubList<char> dataSlice
                    (
                        recvData,
                        sliceOffsets[proci+1]-sliceOffsets[proci],
                        sliceOffsets[proci]
                    );

                    blockOffset[proci] =
                        decomposedBlockData::writeBlockEntry
                        (
                            os,
                            proci,
                            dataSlice
                        );
                }
            }

            startProc += nSendProcs;
        }

        if (UPstream::master(comm))
        {
            ok = osPtr->good();
        }
    }

    if (syncReturnState)
    {
        //- Enable to get synchronised error checking. Is the one that keeps
        //  slaves as slow as the master (which does all the writing)
        Pstream::broadcast(ok, comm);
    }

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<ISstream> isPtr;
    fileName objPath(fileHandler().filePath(false, *this, word::null));
    if (UPstream::master(comm_))
    {
        isPtr.reset(new IFstream(objPath));
        IOobject::readHeader(*isPtr);
    }

    return readBlocks(comm_, isPtr, contentData_, commsType_);
}


bool Foam::decomposedBlockData::writeData(Ostream& os) const
{
    IOobject io(*this);
    IOstreamOption streamOpt(os);

    int verValue;
    int fmtValue;
    // Unneeded: word masterName(name());
    fileName masterLocation(instance()/db().dbDir()/local());

    // Re-read my own data to find out the header information
    if (Pstream::master(comm_))
    {
        ISpanStream headerStream(contentData_);
        io.readHeader(headerStream);

        verValue = headerStream.version().canonical();
        fmtValue = static_cast<int>(headerStream.format());
    }

    // Broadcast header information
    Pstream::broadcasts
    (
        comm_,
        verValue,
        fmtValue,
        // Unneeded: masterName
        io.headerClassName(),
        io.note(),
        // Unneeded: io.instance()
        // Unneeded: io.local()
        masterLocation
    );

    streamOpt.version(IOstreamOption::versionNumber::canonical(verValue));
    streamOpt.format(IOstreamOption::streamFormat(fmtValue));

    if (!Pstream::master(comm_))
    {
        decomposedBlockData::writeHeader
        (
            os,
            streamOpt,  // streamOpt for data
            io.headerClassName(),
            io.note(),
            masterLocation,
            name(),
            dictionary()
        );
    }

    // Write the character data
    if (isA<OFstream>(os))
    {
        // Serial file output - can use writeRaw()
        os.writeRaw(contentData_.cdata(), contentData_.size_bytes());
    }
    else
    {
        // Other cases are less fortunate, and no std::string_view
        std::string str(contentData_.cdata(), contentData_.size_bytes());
        os.writeQuoted(str, false);
    }

    if (!Pstream::master(comm_))
    {
        IOobject::writeEndDivider(os);
    }

    return os.good();
}


bool Foam::decomposedBlockData::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm_))
    {
        // Note: always write binary. These are strings so readable anyway.
        //       They have already be tokenised on the sending side.

        osPtr.reset(new OFstream(objectPath(), IOstreamOption::BINARY));

        // Update meta-data for current state
        const_cast<regIOobject&>
        (
            static_cast<const regIOobject&>(*this)
        ).updateMetaData();

        decomposedBlockData::writeHeader
        (
            *osPtr,
            streamOpt,  // streamOpt for data
            static_cast<const IOobject&>(*this)
        );
    }

    labelList recvSizes;
    gather(comm_, label(contentData_.size_bytes()), recvSizes);

    List<std::streamoff> blockOffsets;
    PtrList<SubList<char>> slaveData;  // dummy slave data
    return writeBlocks
    (
        comm_,
        osPtr,
        blockOffsets,
        contentData_,
        recvSizes,
        slaveData,
        commsType_
    );
}


// ************************************************************************* //
