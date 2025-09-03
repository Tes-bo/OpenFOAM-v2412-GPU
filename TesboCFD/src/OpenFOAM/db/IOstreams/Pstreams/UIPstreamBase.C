/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "error.H"
#include "UIPstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Convert a single character to a word with length 1
inline static Foam::word charToWord(char c)
{
    return Foam::word(std::string(1, c), false);
}


// Adjust stream format based on the flagMask
inline static void processFlags(Istream& is, int flagMask)
{
    if ((flagMask & token::ASCII))
    {
        is.format(IOstreamOption::ASCII);
    }
    else if ((flagMask & token::BINARY))
    {
        is.format(IOstreamOption::BINARY);
    }
}


// Return the position with word boundary alignment
inline static label byteAlign(const label pos, const size_t align)
{
    return
    (
        (align > 1)
      ? (align + ((pos - 1) & ~(align - 1)))
      : pos
    );
}


// Read into compound token (assumed to be a known type)
inline static bool readCompoundToken
(
    token& tok,
    const word& compoundType,
    Istream& is
)
{
    // The isCompound() check is not needed (already checked by caller)
    return tok.readCompoundToken(compoundType, is);
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::UIPstreamBase::checkEof()
{
    if (recvBufPos_ == messageSize_)
    {
        setEof();
    }
}


inline void Foam::UIPstreamBase::prepareBuffer(const size_t align)
{
    recvBufPos_ = byteAlign(recvBufPos_, align);
}


template<class T>
inline void Foam::UIPstreamBase::readFromBuffer(T& val)
{
    prepareBuffer(sizeof(T));

    val = reinterpret_cast<T&>(recvBuf_[recvBufPos_]);
    recvBufPos_ += sizeof(T);
    checkEof();
}


inline void Foam::UIPstreamBase::readFromBuffer
(
    void* data,
    const size_t count
)
{
    if (data)
    {
        const char* const __restrict__ buf = &recvBuf_[recvBufPos_];
        char* const __restrict__ output = reinterpret_cast<char*>(data);

        for (size_t i = 0; i < count; ++i)
        {
            output[i] = buf[i];
        }
    }

    recvBufPos_ += count;
    checkEof();
}


inline Foam::Istream& Foam::UIPstreamBase::readString(std::string& str)
{
    // Use std::string::assign() to copy content, including embedded nul chars.
    // Stripping (when desired) is the responsibility of the sending side.

    size_t len;
    readFromBuffer(len);

    if (len)
    {
        str.assign(&recvBuf_[recvBufPos_], len);
        recvBufPos_ += len;
        checkEof();
    }
    else
    {
        str.clear();
    }

    return *this;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UIPstreamBase::UIPstreamBase
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    DynamicList<char>& receiveBuf,
    label& receiveBufPosition,
    const int tag,
    const label comm,
    const bool clearAtEnd,
    IOstreamOption::streamFormat fmt
)
:
    UPstream(commsType),
    Istream(fmt),
    fromProcNo_(fromProcNo),
    tag_(tag),
    comm_(comm),
    messageSize_(0),
    storedRecvBufPos_(0),
    clearAtEnd_(clearAtEnd),
    recvBuf_(receiveBuf),
    recvBufPos_(receiveBufPosition)
{
    setOpened();
    setGood();
}


Foam::UIPstreamBase::UIPstreamBase
(
    const int fromProcNo,
    PstreamBuffers& buffers
)
:
    UPstream(buffers.commsType()),
    Istream(buffers.format()),
    fromProcNo_(fromProcNo),
    tag_(buffers.tag()),
    comm_(buffers.comm()),
    messageSize_(0),
    storedRecvBufPos_(0),
    clearAtEnd_(buffers.allowClearRecv()),
    recvBuf_(buffers.accessRecvBuffer(fromProcNo)),
    recvBufPos_(buffers.accessRecvPosition(fromProcNo))
{
    if
    (
        commsType() != UPstream::commsTypes::scheduled
     && !buffers.finished()
    )
    {
        FatalErrorInFunction
            << "PstreamBuffers::finishedSends() never called." << endl
            << "Please call PstreamBuffers::finishedSends() after doing"
            << " all your sends (using UOPstream) and before doing any"
            << " receives (using UIPstream)" << Foam::exit(FatalError);
    }

    setOpened();
    setGood();
}


Foam::UIPstreamBase::UIPstreamBase
(
    const DynamicList<char>& receiveBuf,
    IOstreamOption::streamFormat fmt
)
:
    UPstream(UPstream::commsTypes::nonBlocking), // placeholder
    Istream(fmt),
    fromProcNo_(UPstream::masterNo()),      // placeholder
    tag_(UPstream::msgType()),              // placeholder
    comm_(UPstream::commSelf()),            // placeholder
    messageSize_(receiveBuf.size()),        // Message == buffer
    storedRecvBufPos_(0),
    clearAtEnd_(false),   // Do not clear recvBuf if at end!!
    recvBuf_
    (
        // The receive buffer is never modified with this code path
        const_cast<DynamicList<char>&>(receiveBuf)
    ),
    recvBufPos_(storedRecvBufPos_)          // Internal reference
{
    setOpened();
    setGood();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UIPstreamBase::~UIPstreamBase()
{
    if (clearAtEnd_ && eof())
    {
        if (debug)
        {
            Perr<< "UIPstreamBase Destructor : tag:" << tag_
                << " fromProcNo:" << fromProcNo_
                << " clearing receive buffer of size "
                << recvBuf_.size()
                << " messageSize_:" << messageSize_ << endl;
        }
        recvBuf_.clearStorage();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Istream& Foam::UIPstreamBase::read(token& t)
{
    // Return the put back token if it exists
    // - with additional handling for special stream flags
    if (Istream::getBack(t))
    {
        if (t.isFlag())
        {
            processFlags(*this, t.flagToken());
        }
        else
        {
            return *this;
        }
    }


    // Reset token, adjust its line number according to the stream
    t.reset();
    t.lineNumber(this->lineNumber());


    // Read character, return on error
    // - with additional handling for special stream flags

    char c;
    do
    {
        if (!read(c))
        {
            t.setBad();   // Error
            return *this;
        }

        if (c == token::FLAG)
        {
            char flagVal;

            if (read(flagVal))
            {
                processFlags(*this, flagVal);
            }
            else
            {
                t.setBad();   // Error
                return *this;
            }
        }
    }
    while (c == token::FLAG);


    // Analyse input starting with this character.
    switch (c)
    {
        // Punctuation
        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::PLUS :
        case token::MINUS :
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // The word-variants
        case token::tokenType::WORD :
        case token::tokenType::DIRECTIVE :
        {
            word val;
            if (readString(val))
            {
                if
                (
                    !token::compound::isCompound(val)
                 || !readCompoundToken(t, val, *this)
                )
                {
                    t = std::move(val);
                    t.setType(token::tokenType(c));
                }
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // The string-variants
        case token::tokenType::STRING :
        case token::tokenType::EXPRESSION :
        case token::tokenType::VARIABLE :
        case token::tokenType::VERBATIM :
        case token::tokenType::CHAR_DATA :
        {
            string val;
            if (readString(val))
            {
                t = std::move(val);
                t.setType(token::tokenType(c));
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Label
        case token::tokenType::LABEL :
        {
            label val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Float
        case token::tokenType::FLOAT :
        {
            float val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Double
        case token::tokenType::DOUBLE :
        {
            double val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Character (returned as a single character word) or error
        default:
        {
            if (isalpha(c))
            {
                t = charToWord(c);
                return *this;
            }

            setBad();
            t.setBad();

            return *this;
        }
    }
}


Foam::Istream& Foam::UIPstreamBase::read(char& c)
{
    c = recvBuf_[recvBufPos_];
    ++recvBufPos_;
    checkEof();
    return *this;
}


Foam::Istream& Foam::UIPstreamBase::read(word& str)
{
    return readString(str);
}


Foam::Istream& Foam::UIPstreamBase::read(string& str)
{
    return readString(str);
}


Foam::Istream& Foam::UIPstreamBase::read(label& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstreamBase::read(float& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstreamBase::read(double& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstreamBase::read(char* data, std::streamsize count)
{
    if (count)
    {
        // For count == 0, a no-op
        // - see UOPstream::write(const char*, streamsize)
        beginRawRead();
        readRaw(data, count);
        endRawRead();
    }

    return *this;
}


Foam::Istream& Foam::UIPstreamBase::readRaw(char* data, std::streamsize count)
{
    // No check for IOstreamOption::BINARY since this is either done in the
    // beginRawRead() method, or the caller knows what they are doing.

    // Any alignment must have been done prior to this call
    readFromBuffer(data, count);
    return *this;
}


bool Foam::UIPstreamBase::beginRawRead()
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    // Align on word boundary (64-bit)
    // - as per read(const char*, streamsize)
    // The check for zero-size will have been done by the caller
    prepareBuffer(8);

    return true;
}


// Not needed yet
///
/// //- The current get position (tellg) in the buffer
/// label pos() const;
///
/// Foam::label Foam::UIPstreamBase::pos() const
/// {
///     return recvBufPos_;
/// }


void Foam::UIPstreamBase::rewind()
{
    recvBufPos_ = 0;  // Assume the entire buffer is for us to read from
    setOpened();
    setGood();
    if (recvBuf_.empty() || !messageSize_)
    {
        setEof();
    }
}


void Foam::UIPstreamBase::print(Ostream& os) const
{
    os  << "Reading from processor " << fromProcNo_
        << " using communicator " << comm_
        <<  " and tag " << tag_ << Foam::endl;
}


// ************************************************************************* //
