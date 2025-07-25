/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::UIPstreamBase

Description
    Base class for input inter-processor communications stream
    (ie, parallel streams).
    Not to be used directly, thus contructors are protected.

SourceFiles
    UIPstreamBase.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef Foam_UIPstream_H
#define Foam_UIPstream_H

#include "UPstream.H"
#include "Istream.H"
#include "DynamicList.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class UIPstreamBase Declaration
\*---------------------------------------------------------------------------*/

class UIPstreamBase
:
    public UPstream,
    public Istream
{
    // Private Member Functions

        //- Check buffer position against messageSize_ for EOF
        inline void checkEof();

        //- Prepare receive buffer by adjusting alignment
        inline void prepareBuffer(const size_t align);

        //- Read a T from the receive buffer
        template<class T>
        inline void readFromBuffer(T& val);

        //- Read count bytes of data from the receive buffer.
        //  Prior data alignment is done by prepareBuffer
        //  Reading into a null pointer behaves like a forward seek
        inline void readFromBuffer(void* data, const size_t count);

        //- Read string length and string content
        inline Istream& readString(std::string& str);


protected:

    // Protected Data

        //- Source rank for the data
        const int fromProcNo_;

        //- Message tag for communication
        const int tag_;

        //- The communicator index
        const int comm_;

        //- The message size, read on bufferIPCrecv or set directly
        label messageSize_;

        //- Receive position in buffer data, if ony
        //- If there is no external location for recvBufPos_
        label storedRecvBufPos_;

        //- Clear the receive buffer on termination (in the destructor)
        const bool clearAtEnd_;

        //- Reference to the receive buffer data
        DynamicList<char>& recvBuf_;

        //- Reference to the receive position in buffer data
        label& recvBufPos_;


    // Protected Constructors

        //- Construct given process index to read from using the given
        //- attached receive buffer, optional communication characteristics
        //- and IO format
        UIPstreamBase
        (
            const UPstream::commsTypes commsType,
            const int fromProcNo,
            DynamicList<char>& receiveBuf,
            label& receiveBufPosition,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool clearAtEnd = false,   // destroy receiveBuf if at end
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct given buffers
        UIPstreamBase(const int fromProcNo, PstreamBuffers& buffers);

        //- Construct for an externally obtained buffer.
        //  The parameter is allowed to be const (since reading will not
        //  affect it), but must reference a concrete variable.
        UIPstreamBase
        (
            const DynamicList<char>& receiveBuf,
            IOstreamOption::streamFormat fmt
        );

public:

    //- Destructor. Optionally clears external receive buffer.
    virtual ~UIPstreamBase();


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


    // Read Functions

        //- Return next token from stream
        virtual Istream& read(token&) override;

        //- Read a character
        virtual Istream& read(char& c) override;

        //- Read a word
        virtual Istream& read(word& str) override;

        // Read a string
        virtual Istream& read(string& str) override;

        //- Read a label
        virtual Istream& read(label& val) override;

        //- Read a float
        virtual Istream& read(float& val) override;

        //- Read a double
        virtual Istream& read(double& val) override;

        //- Read binary block with 8-byte alignment.
        //- Reading into a null pointer behaves like a forward seek of
        //- count characters.
        virtual Istream& read(char* data, std::streamsize count) override;

        //- Low-level raw binary read.
        //- Reading into a null pointer behaves like a forward seek of
        //- count characters.
        virtual Istream& readRaw(char* data, std::streamsize count) override;

        //- Start of low-level raw binary read
        virtual bool beginRawRead() override;

        //- End of low-level raw binary read
        virtual bool endRawRead() override { return true; }


    // Positioning

        //- Rewind the receive stream position so that it may be read again
        virtual void rewind() override;


    // Print

        //- Print stream description to Ostream
        void print(Ostream& os) const override;
};


/*---------------------------------------------------------------------------*\
                           Class UIPstream Declaration
\*---------------------------------------------------------------------------*/

//- Input inter-processor communications stream
//- using MPI send/recv etc. - operating on external buffer.
class UIPstream
:
    public UIPstreamBase
{
    // Private Member Functions

        //- Initial buffer recv, called by constructor (blocking | scheduled)
        void bufferIPCrecv();


public:

    // Constructors

        //- Construct given process index to read from using the given
        //- attached receive buffer, optional communication characteristics
        //- and IO format
        UIPstream
        (
            const UPstream::commsTypes commsType,
            const int fromProcNo,
            DynamicList<char>& receiveBuf,
            label& receiveBufPosition,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool clearAtEnd = false,   // destroy receiveBuf if at end
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct given buffers
        UIPstream(const int fromProcNo, PstreamBuffers& buffers);

        //- Construct for reading from a standalone buffer that has
        //- been obtained externally by the caller.
        //  The parameter is allowed to be const (since reading will not
        //  affect it), but must reference a concrete variable.
        explicit UIPstream
        (
            const DynamicList<char>& recvBuf,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    //- Destructor
    virtual ~UIPstream() = default;


    // Member Functions

        //- Use all read methods from base
        using UIPstreamBase::read;


    // Static Functions

        //- Read buffer contents from given processor.
        //  \return the message size (bytes read). May change in the future
        static std::streamsize read
        (
            const UPstream::commsTypes commsType,
            const int fromProcNo,
            char* buf,
            const std::streamsize bufSize,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        );

        //- Read buffer contents (non-blocking) from given processor.
        //  \return the message size (bytes read). May change in the future
        inline static std::streamsize read
        (
            //! [out] request information
            UPstream::Request& req,
            const int fromProcNo,
            char* buf,
            const std::streamsize bufSize,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        )
        {
            return UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                fromProcNo,
                buf,
                bufSize,
                tag,
                comm,
               &req
            );
        }

        //- Read into UList storage from given processor.
        //  Only valid for contiguous data types.
        //  \return the message size (bytes read). May change in the future
        template<class Type>
        inline static std::streamsize read
        (
            const UPstream::commsTypes commsType,
            const int fromProcNo,
            UList<Type>& buffer,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        )
        {
            return UIPstream::read
            (
                commsType,
                fromProcNo,
                buffer.data_bytes(),
                buffer.size_bytes(),
                tag,
                comm,
                req
            );
        }

        //- Read into SubList storage from given processor.
        //  Only valid for contiguous data types.
        //  \return the message size (bytes read). May change in the future
        template<class Type>
        inline static std::streamsize read
        (
            const UPstream::commsTypes commsType,
            const int fromProcNo,
            SubList<Type> buffer,  // passed by shallow copy
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        )
        {
            return UIPstream::read
            (
                commsType,
                fromProcNo,
                buffer.data_bytes(),
                buffer.size_bytes(),
                tag,
                comm,
                req
            );
        }

        //- Read into UList storage (non-blocking) from given processor.
        //  Only valid for contiguous data types.
        //  \return the message size (bytes read). May change in the future
        template<class Type>
        inline static std::streamsize read
        (
            //! [out] request information
            UPstream::Request& req,
            const int fromProcNo,
            UList<Type>& buffer,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        )
        {
            return UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                fromProcNo,
                buffer.data_bytes(),
                buffer.size_bytes(),
                tag,
                comm,
               &req
            );
        }

        //- Read into SubList storage (non-blocking) from given processor.
        //  Only valid for contiguous data types.
        //  \return the message size (bytes read). May change in the future
        template<class Type>
        inline static std::streamsize read
        (
            //! [out] request information
            UPstream::Request& req,
            const int fromProcNo,
            SubList<Type> buffer,  // passed by shallow copy
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        )
        {
            return UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                fromProcNo,
                buffer.data_bytes(),
                buffer.size_bytes(),
                tag,
                comm,
               &req
            );
        }
};


/*---------------------------------------------------------------------------*\
                         Class UIPBstream Declaration
\*---------------------------------------------------------------------------*/

//- Input inter-processor communications stream
//- using MPI broadcast - operating on external buffer.
class UIPBstream
:
    public UIPstreamBase
{
    // Private Member Functions

        //- Initial buffer recv via broadcast, called by constructor
        void bufferIPCrecv();


public:

    // Constructors

        //- Construct given process index to read from using the given
        //- attached receive buffer, optional communication characteristics
        //- and IO format
        UIPBstream
        (
            const UPstream::commsTypes,     //!< irrelevant
            const int rootProcNo,           //!< normally UPstream::masterNo()
            DynamicList<char>& receiveBuf,
            label& receiveBufPosition,
            const int tag = UPstream::msgType(),  //!< irrelevant
            const label comm = UPstream::worldComm,
            const bool clearAtEnd = false,  //!< destroy receiveBuf if at end
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    //- Destructor
    virtual ~UIPBstream() = default;


    // Member Functions

        //- Use all read methods from base
        using UIPstreamBase::read;


    // Static Functions

        //- Wrapped version of UPstream::broadcast
        //  \return the message size (bytes read). May change in the future
        static std::streamsize read
        (
            const int rootProcNo,  //!< normally UPstream::masterNo()
            char* buf,
            const std::streamsize bufSize,
            const label comm = UPstream::worldComm
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
