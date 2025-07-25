/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::UOPstreamBase

Description
    Base class for output inter-processor communications stream
    (ie, parallel streams).
    Not to be used directly, thus contructors are protected.

SourceFiles
    UOPstreamBase.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef Foam_UOPstream_H
#define Foam_UOPstream_H

#include "UPstream.H"
#include "Ostream.H"
#include "DynamicList.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class UOPstreamBase Declaration
\*---------------------------------------------------------------------------*/

class UOPstreamBase
:
    public UPstream,
    public Ostream
{
    // Private Member Functions

        //- Prepare send buffer for count bytes of output,
        //- with specified alignment.
        inline void prepareBuffer(const size_t count, const size_t align);

        //- Generic write data to the send buffer, aligned by sizeof(T)
        template<class T>
        inline void writeToBuffer(const T& val);

        //- Write count bytes of data to the send buffer,
        //- using align byte alignment
        inline void writeToBuffer
        (
            const void* data,
            const size_t count,
            const size_t align
        );

        //- Add a single char to the send buffer. No alignment needed
        inline void putChar(const char c);

        //- Write string length and content, including embedded nul chars.
        inline void putString(const char* str, const size_t len);

        //- Write string length and content, including embedded nul chars.
        inline void putString(const std::string& str);


protected:

    // Protected Data

        //- Destination rank for the data
        const int toProcNo_;

        //- Message tag for communication
        const int tag_;

        //- The communicator index
        const int comm_;

        //- Call bufferIPCsend on termination (in the destructor)
        bool sendAtDestruct_;

        //- Reference to the send buffer data
        DynamicList<char>& sendBuf_;


    // Protected Constructors

        //- Construct given process index to write to using the given
        //- attached send buffer, optional communication characteristics
        //- and IO format
        UOPstreamBase
        (
            const UPstream::commsTypes commsType,
            const int toProcNo,
            DynamicList<char>& sendBuf,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool sendAtDestruct = true,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct given buffers
        UOPstreamBase(const int toProcNo, PstreamBuffers& buffers);

        //- Construct for externally obtained buffers
        UOPstreamBase
        (
            DynamicList<char>& sendBuf,
            IOstreamOption::streamFormat fmt
        );

public:

    //- Destructor.
    virtual ~UOPstreamBase();


    // Member Functions

    // Stream State Functions

        //- Return current stream flags.
        //- Dummy for parallel stream, returns 0.
        virtual std::ios_base::fmtflags flags() const override
        {
            return std::ios_base::fmtflags(0);
        }

        //- Set stream flags, return old stream flags.
        //- Dummy for parallel stream, returns 0.
        virtual std::ios_base::fmtflags flags(std::ios_base::fmtflags) override
        {
            return std::ios_base::fmtflags(0);
        }


    // Write Functions

        //- Inherit write methods from Ostream
        using Ostream::writeQuoted;

        //- Write token to stream or otherwise handle it.
        //  \return false if the token type was not handled by this method
        virtual bool write(const token& tok) override;

        //- Write single character. Whitespace is suppressed.
        virtual Ostream& write(const char c) override;

        //- Write character/string content, with/without surrounding quotes
        virtual Ostream& writeQuoted
        (
            const char* str,
            std::streamsize len,
            const bool quoted=true
        ) override;

        //- Write the word-characters of a character string.
        //  Sends as a single char, or as word.
        virtual Ostream& write(const char* str) override;

        //- Write word
        virtual Ostream& write(const word& str) override;

        //- Write string
        virtual Ostream& write(const std::string& str) override;

        //- Write int32_t as a label
        virtual Ostream& write(const int32_t val) override;

        //- Write int64_t as a label
        virtual Ostream& write(const int64_t val) override;

        //- Write float
        virtual Ostream& write(const float val) override;

        //- Write double
        virtual Ostream& write(const double val) override;

        //- Write binary block with 8-byte alignment.
        virtual Ostream& write
        (
            const char* data,
            std::streamsize count
        ) override;

        //- Low-level raw binary output.
        virtual Ostream& writeRaw
        (
            const char* data,
            std::streamsize count
        ) override;

        //- Begin marker for low-level raw binary output.
        //  The count indicates the number of bytes for subsequent
        //  writeRaw calls.
        virtual bool beginRawWrite(std::streamsize count) override;

        //- End marker for low-level raw binary output.
        virtual bool endRawWrite() override
        {
            return true;
        }

        //- Add indentation characters
        virtual void indent() override
        {}


    // Stream state functions

        //- Flush stream
        virtual void flush() override
        {}

        //- Add newline and flush stream
        virtual void endl() override
        {}

        //- Get the current padding character
        //  \return previous padding character
        virtual char fill() const override
        {
            return 0;
        }

        //- Set padding character for formatted field up to field width
        virtual char fill(const char) override
        {
            return 0;
        }

        //- Get width of output field
        virtual int width() const override
        {
            return 0;
        }

        //- Set width of output field
        //  \return previous width
        virtual int width(const int) override
        {
             return 0;
        }

        //- Get precision of output field
        virtual int precision() const override
        {
             return 0;
        }

        //- Set precision of output field
        //  \return old precision
        virtual int precision(const int) override
        {
             return 0;
        }


    // Positioning

        //- Rewind the send buffer for overwriting
        virtual void rewind();


    // Print

        //- Print stream description to Ostream
        void print(Ostream& os) const override;
};


/*---------------------------------------------------------------------------*\
                          Class UOPstream Declaration
\*---------------------------------------------------------------------------*/

//- Output inter-processor communications stream
//- using MPI send/recv etc. - operating on external buffer.
class UOPstream
:
    public UOPstreamBase
{
    // Private Member Functions

        //- Buffer send, usually called by destructor
        bool bufferIPCsend();


public:

    // Constructors

        //- Construct given process index to write to using the given
        //- attached send buffer, optional communication characteristics
        //- and IO format
        UOPstream
        (
            const UPstream::commsTypes commsType,
            const int toProcNo,
            DynamicList<char>& sendBuf,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool sendAtDestruct = true,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct given buffers
        UOPstream(const int toProcNo, PstreamBuffers& buffers);

        //- Construct for writing into a standalone buffer.
        //- Data transfer is handled externally by the caller.
        explicit UOPstream
        (
            DynamicList<char>& sendBuf,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    //- Destructor, usually sends buffer on destruct.
    virtual ~UOPstream();


    // Member Functions

        //- Send buffer contents now and not in destructor [advanced usage].
        //- Returns true on success
        bool send();

        //- Use all write methods from base class
        using UOPstreamBase::write;


    // Static Functions

        //- Write buffer contents to given processor
        //  \return True on success
        static bool write
        (
            const UPstream::commsTypes commsType,
            const int toProcNo,
            const char* buf,
            const std::streamsize bufSize,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr,
            const UPstream::sendModes sendMode = UPstream::sendModes::normal
        );

        //- Write buffer contents (non-blocking) to given processor
        //  \return True on success
        inline static bool write
        (
            //! [out] request information
            UPstream::Request& req,
            const int toProcNo,
            const char* buf,
            const std::streamsize bufSize,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const UPstream::sendModes sendMode = UPstream::sendModes::normal
        )
        {
            return UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                toProcNo,
                buf,
                bufSize,
                tag,
                comm,
               &req,
                sendMode
            );
        }

        //- Write UList contents to given processor.
        //  Only valid for contiguous data types.
        //  \return True on success
        template<class Type>
        inline static bool write
        (
            const UPstream::commsTypes commsType,
            const int toProcNo,
            const UList<Type>& buffer,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr,
            const UPstream::sendModes sendMode = UPstream::sendModes::normal
        )
        {
            return UOPstream::write
            (
                commsType,
                toProcNo,
                buffer.cdata_bytes(),
                buffer.size_bytes(),
                tag,
                comm,
                req,
                sendMode
            );
        }

        //- Write UList contents (non-blocking) to given processor.
        //  Only valid for contiguous data types.
        //  \return True on success
        template<class Type>
        inline static bool write
        (
            //! [out] request information
            UPstream::Request& req,
            const int toProcNo,
            const UList<Type>& buffer,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const UPstream::sendModes sendMode = UPstream::sendModes::normal
        )
        {
            return UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                toProcNo,
                buffer.cdata_bytes(),
                buffer.size_bytes(),
                tag,
                comm,
               &req,
                sendMode
            );
        }
};


/*---------------------------------------------------------------------------*\
                         Class UOPBstream Declaration
\*---------------------------------------------------------------------------*/

//- Output inter-processor communications stream
//- using MPI broadcast - operating on external buffer.
//
// \note does not use commsType, tag etc.
class UOPBstream
:
    public UOPstreamBase
{
    // Private Member Functions

        //- Buffer send, usually called by destructor
        bool bufferIPCsend();


public:

    // Constructors

        //- Construct given process index to write to using the given
        //- attached send buffer, optional communication characteristics
        //- and IO format
        UOPBstream
        (
            const UPstream::commsTypes,     //!< irrelevant
            const int toProcNo,             //!< normally UPstream::masterNo()
            DynamicList<char>& sendBuf,
            const int tag = UPstream::msgType(),  //!< irrelevant
            const label comm = UPstream::worldComm,
            const bool sendAtDestruct = true,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    //- Destructor, usually sends buffer on destruct.
    virtual ~UOPBstream();


    // Member Functions

        //- Use all write methods from base
        using UOPstreamBase::write;


    // Static Functions

        //- Wrapped version of UPstream::broadcast with const-cast
        //  \return True on success
        static bool write
        (
            const int rootProcNo,  //!< normally UPstream::masterNo()
            const char* buf,
            const std::streamsize bufSize,
            const label comm = UPstream::worldComm
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
