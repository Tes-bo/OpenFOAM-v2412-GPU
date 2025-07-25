/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "UOPstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

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

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::UOPstreamBase::prepareBuffer
(
    const size_t count,
    const size_t align
)
{
    if (!count)
    {
        return;
    }

    // Align for the next output position
    const label pos = byteAlign(sendBuf_.size(), align);

    // Extend buffer (as required)
    sendBuf_.reserve(Foam::max(label(1024), label(pos + count)));

    // Move to the aligned output position. Fill any gap with nul char.
    sendBuf_.resize(pos, '\0');
}


template<class T>
inline void Foam::UOPstreamBase::writeToBuffer(const T& val)
{
    writeToBuffer(&val, sizeof(T), sizeof(T));
}


inline void Foam::UOPstreamBase::writeToBuffer
(
    const void* data,
    const size_t count,
    const size_t align
)
{
    if (!count)
    {
        return;
    }

    prepareBuffer(count, align);

    // The aligned output position
    const label pos = sendBuf_.size();

    // Extend the addressable range for direct pointer access
    sendBuf_.resize(pos + count);

    char* const __restrict__ buf = (sendBuf_.data() + pos);
    const char* const __restrict__ input = reinterpret_cast<const char*>(data);

    for (size_t i = 0; i < count; ++i)
    {
        buf[i] = input[i];
    }
}


inline void Foam::UOPstreamBase::putChar(const char c)
{
    if (!sendBuf_.capacity())
    {
        sendBuf_.setCapacity(1024);
    }
    sendBuf_.push_back(c);
}


inline void Foam::UOPstreamBase::putString
(
    const char* str,
    const size_t len
)
{
    writeToBuffer(len);
    writeToBuffer(str, len, 1);  // no-op when len == 0
}


inline void Foam::UOPstreamBase::putString(const std::string& str)
{
    putString(str.data(), str.size());
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UOPstreamBase::UOPstreamBase
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    DynamicList<char>& sendBuf,
    const int tag,
    const label comm,
    const bool sendAtDestruct,
    IOstreamOption::streamFormat fmt
)
:
    UPstream(commsType),
    Ostream(fmt),
    toProcNo_(toProcNo),
    tag_(tag),
    comm_(comm),
    sendAtDestruct_(sendAtDestruct),
    sendBuf_(sendBuf)
{
    setOpened();
    setGood();
}


Foam::UOPstreamBase::UOPstreamBase(const int toProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType()),
    Ostream(buffers.format()),
    toProcNo_(toProcNo),
    tag_(buffers.tag()),
    comm_(buffers.comm()),
    sendAtDestruct_(buffers.commsType() != UPstream::commsTypes::nonBlocking),
    sendBuf_(buffers.accessSendBuffer(toProcNo))
{
    setOpened();
    setGood();
}


Foam::UOPstreamBase::UOPstreamBase
(
    DynamicList<char>& sendBuf,
    IOstreamOption::streamFormat fmt
)
:
    UPstream(UPstream::commsTypes::nonBlocking), // placeholder
    Ostream(fmt),
    toProcNo_(UPstream::masterNo()),        // placeholder
    tag_(UPstream::msgType()),              // placeholder
    comm_(UPstream::commSelf()),            // placeholder
    sendAtDestruct_(false),   // Never sendAtDestruct!!
    sendBuf_(sendBuf)
{
    sendBuf_.clear();  // Overwrite into buffer
    setOpened();
    setGood();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UOPstreamBase::~UOPstreamBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::UOPstreamBase::write(const token& tok)
{
    // Direct token handling only for some types

    switch (tok.type())
    {
        case token::tokenType::FLAG :
        {
            putChar(token::tokenType::FLAG);
            putChar(tok.flagToken());

            return true;
        }

        // The word-variants
        case token::tokenType::WORD :
        case token::tokenType::DIRECTIVE :
        {
            putChar(tok.type());
            putString(tok.wordToken());

            return true;
        }

        // The string-variants
        case token::tokenType::STRING :
        case token::tokenType::EXPRESSION :
        case token::tokenType::VARIABLE :
        case token::tokenType::VERBATIM :
        case token::tokenType::CHAR_DATA :
        {
            putChar(tok.type());
            putString(tok.stringToken());

            return true;
        }

        default:
            break;
    }

    return false;
}


Foam::Ostream& Foam::UOPstreamBase::write(const char c)
{
    if (!isspace(c))
    {
        putChar(c);
    }

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::writeQuoted
(
    const char* str,
    std::streamsize len,
    const bool quoted
)
{
    if (quoted)
    {
        putChar(token::tokenType::STRING);
    }
    else
    {
        putChar(token::tokenType::WORD);
    }
    putString(str, len);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const char* str)
{
    const word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        return write(nonWhiteChars[0]);
    }
    else if (nonWhiteChars.size())
    {
        return write(nonWhiteChars);
    }

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const word& str)
{
    putChar(token::tokenType::WORD);
    putString(str);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const std::string& str)
{
    putChar(token::tokenType::STRING);
    putString(str);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const int32_t val)
{
    putChar(token::tokenType::LABEL);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const int64_t val)
{
    putChar(token::tokenType::LABEL);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const float val)
{
    putChar(token::tokenType::FLOAT);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const double val)
{
    putChar(token::tokenType::DOUBLE);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write
(
    const char* data,
    std::streamsize count
)
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    // Align on word boundary (64-bit)
    writeToBuffer(data, count, 8);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::writeRaw
(
    const char* data,
    std::streamsize count
)
{
    // No check for IOstreamOption::BINARY since this is either done in the
    // beginRawWrite() method, or the caller knows what they are doing.

    // Previously aligned and sizes reserved via beginRawWrite()
    writeToBuffer(data, count, 1);

    return *this;
}


bool Foam::UOPstreamBase::beginRawWrite(std::streamsize count)
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    // Align on word boundary (64-bit)
    // - as per write(const char*, streamsize)
    prepareBuffer(count, 8);

    return true;
}


// Not needed yet
///
/// //- The current put position (tellp) in the buffer
/// label pos() const;
///
/// Foam::label Foam::UOPstreamBase::pos() const
/// {
///     return sendBuf_.size();
/// }


void Foam::UOPstreamBase::rewind()
{
    sendBuf_.clear();  // Overwrite into buffer
    setOpened();
    setGood();
}


void Foam::UOPstreamBase::print(Ostream& os) const
{
    os  << "Writing to processor " << toProcNo_
        << " from processor " << myProcNo() << " in communicator " << comm_
        << " and tag " << tag_ << Foam::endl;
}


// ************************************************************************* //
