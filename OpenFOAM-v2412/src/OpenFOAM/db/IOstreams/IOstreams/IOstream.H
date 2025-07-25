/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::IOstream

Description
    An IOstream is an abstract base class for all input/output systems; be
    they streams, files, token lists etc.

    The basic operations are construct, close, read token, read primitive
    and read binary block.  In addition version control and line number
    counting is incorporated.  Usually one would use the read primitive
    member functions, but if one were reading a stream on unknown data
    sequence one can read token by token, and then analyse.

SourceFiles
    IOstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOstream_H
#define Foam_IOstream_H

#include "char.H"
#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "scalar.H"
#include "fileName.H"
#include "InfoProxy.H"
#include "IOstreamOption.H"

#include <iostream>

using std::ios_base;
using std::istream;
using std::ostream;

using std::cerr;
using std::cin;
using std::cout;

// Additional constructors and methods (as per v2012 and earlier)
#define Foam_IOstream_extras
// COMPAT_OPENFOAM_ORG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOstream;

template<>
Ostream& operator<<(Ostream&, const InfoProxy<IOstream>&);


/*---------------------------------------------------------------------------*\
                           Class IOstream Declaration
\*---------------------------------------------------------------------------*/

class IOstream
:
    public IOstreamOption
{
public:

    // Public Data Types

        //- Enumeration for stream open/closed state
        enum streamAccess : char
        {
            CLOSED = 0,         //!< The stream is not open
            OPENED              //!< The stream is open
        };


    // Public Static Data

        //- Default precision
        static unsigned int precision_;


protected:

    // Protected Data

        //- Name for any generic stream - normally treat as readonly
        static fileName staticName_;

        //- Mirror of internal stream io state
        std::ios_base::iostate ioState_;

        //- The stream open/closed state
        streamAccess openClosed_;

        //- The sizeof (label), possibly read from the header
        unsigned char sizeofLabel_;

        //- The sizeof (scalar), possibly read from the header
        unsigned char sizeofScalar_;

        //- The file line
        label lineNumber_;


    // Protected Member Functions

    // Access

        //- Set stream opened
        void setOpened() noexcept
        {
            openClosed_ = OPENED;
        }

        //- Set stream closed
        void setClosed() noexcept
        {
            openClosed_ = CLOSED;
        }

        //- Set stream state
        void setState(std::ios_base::iostate state) noexcept
        {
            ioState_ = state;
        }

        //- Set stream state to be good
        void setGood() noexcept
        {
            ioState_ = std::ios_base::goodbit;
        }


public:

    // Generated Methods

        //- Copy construct
        IOstream(const IOstream&) = default;

        //- Destructor
        virtual ~IOstream() = default;


    // Constructors

        //- Default construct (ASCII, uncompressed),
        //- construct with specified stream option
        explicit IOstream(IOstreamOption streamOpt = IOstreamOption())
        :
            IOstreamOption(streamOpt),
            ioState_(std::ios_base::goodbit),
            openClosed_(CLOSED),
            sizeofLabel_(static_cast<unsigned char>(sizeof(label))),
            sizeofScalar_(static_cast<unsigned char>(sizeof(scalar))),
            lineNumber_(0)
        {
            setBad();
        }

        //- Construct with format, version (compression)
        IOstream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            IOstream(IOstreamOption(fmt, ver, cmp))
        {}


    // Member Functions

    // Access

        //- The name of the stream.
        virtual const fileName& name() const;

        //- Return the name of the stream relative to the current case.
        //  Uses argList::envRelativePath()
        fileName relativeName() const;


    // Check

        //- Check IOstream status for given operation.
        //  Print IOstream state or generate a FatalIOError
        //  when an error has occurred.
        //  The base implementation is a fatalCheck
        virtual bool check(const char* operation) const;

        //- Check IOstream status for given operation.
        //  Generate a FatalIOError when an error has occurred.
        bool fatalCheck(const char* operation) const;

        //- True if stream has been opened
        bool opened() const noexcept
        {
            return openClosed_ == OPENED;
        }

        //- True if stream is closed
        bool closed() const noexcept
        {
            return openClosed_ == CLOSED;
        }

        //- True if next operation might succeed
        bool good() const noexcept
        {
            return ioState_ == 0;  // == goodbit
        }

        //- True if end of input seen
        bool eof() const noexcept
        {
            return ioState_ & std::ios_base::eofbit;
        }

        //- True if next operation will fail
        bool fail() const noexcept
        {
            return ioState_ & (std::ios_base::badbit | std::ios_base::failbit);
        }

        //- True if stream is corrupted
        bool bad() const noexcept
        {
            return ioState_ & std::ios_base::badbit;
        }

        //- Return true if the stream has not failed
        explicit operator bool() const noexcept
        {
            return !fail();
        }

        //- Return true if the stream has failed
        bool operator!() const noexcept
        {
            return fail();
        }


    // Element sizes (precision)

        //- The sizeof (label) in bytes associated with the stream
        unsigned labelByteSize() const noexcept
        {
            return static_cast<unsigned>(sizeofLabel_);
        }

        //- The sizeof (scalar) in bytes associated with the stream
        unsigned scalarByteSize() const noexcept
        {
            return static_cast<unsigned>(sizeofScalar_);
        }

        //- Set the sizeof (label) in bytes associated with the stream
        void setLabelByteSize(unsigned nbytes) noexcept
        {
            sizeofLabel_ = static_cast<unsigned char>(nbytes);
        }

        //- Set the sizeof (scalar) in bytes associated with the stream
        void setScalarByteSize(unsigned nbytes) noexcept
        {
            sizeofScalar_ = static_cast<unsigned char>(nbytes);
        }


        //- Check if the label byte-size associated with the stream
        //- is the same as the given type
        template<class T = label>
        typename std::enable_if<std::is_integral<T>::value, bool>::type
        checkLabelSize() const noexcept
        {
            return sizeofLabel_ == sizeof(T);
        }

        //- Check if the scalar byte-size associated with the stream
        //- is the same as the given type
        template<class T = scalar>
        typename std::enable_if<std::is_floating_point<T>::value, bool>::type
        checkScalarSize() const noexcept
        {
            return sizeofScalar_ == sizeof(T);
        }


    // Stream State Functions

        //- Const access to the current stream line number
        label lineNumber() const noexcept
        {
            return lineNumber_;
        }

        //- Non-const access to the current stream line number
        label& lineNumber() noexcept
        {
            return lineNumber_;
        }

        //- Set the stream line number
        //  \return the previous value
        label lineNumber(const label num) noexcept
        {
            label old(lineNumber_);
            lineNumber_ = num;
            return old;
        }

        //- Return the default precision
        static unsigned int defaultPrecision() noexcept
        {
            return precision_;
        }

        //- Reset the default precision
        //  \return the previous default precision
        static unsigned int defaultPrecision(unsigned int prec) noexcept
        {
            unsigned int old(precision_);
            precision_ = prec;
            return old;
        }

        //- Set the minimum default precision
        //  \return the previous default precision
        static unsigned int minPrecision(unsigned int prec) noexcept
        {
            unsigned int old(precision_);
            if (precision_ < prec)
            {
                precision_ = prec;
            }
            return old;
        }

        //- Set stream state as reached 'eof'
        void setEof() noexcept
        {
            ioState_ |= std::ios_base::eofbit;
        }

        //- Set stream state as 'failed'
        void setFail() noexcept
        {
            ioState_ |= std::ios_base::failbit;
        }

        //- Set stream state to be 'bad'
        void setBad() noexcept
        {
            ioState_ |= std::ios_base::badbit;
        }

        //- Return current stream flags
        virtual std::ios_base::fmtflags flags() const = 0;

        //- Set stream flags, return old stream flags
        virtual std::ios_base::fmtflags flags(std::ios_base::fmtflags) = 0;

        //- Set stream flag(s), return old stream flags
        std::ios_base::fmtflags setf(std::ios_base::fmtflags f)
        {
            return flags(flags() | f);
        }

        //- Set stream flag(s) with mask, return old stream flags
        std::ios_base::fmtflags setf
        (
            const std::ios_base::fmtflags f,
            const std::ios_base::fmtflags mask
        )
        {
            return flags((flags() & ~mask) | (f & mask));
        }

        //- Unset flags of stream
        void unsetf(std::ios_base::fmtflags f)
        {
            flags(flags() & ~f);
        }


    // Print

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const;

        //- Print information about the stream state bits
        void print(Ostream& os, const int streamState) const;


    // Info

        //- Return info proxy,
        //- used to print IOstream information to a stream
        InfoProxy<IOstream> info() const noexcept { return *this; }
};


// --------------------------------------------------------------------
// ------ Manipulators (not taking arguments)
// --------------------------------------------------------------------

//- An IOstream manipulator
typedef IOstream& (*IOstreamManip)(IOstream&);

//- operator<< handling for manipulators without arguments
inline IOstream& operator<<(IOstream& io, IOstreamManip f)
{
    return f(io);
}


inline IOstream& dec(IOstream& io)
{
    io.setf(std::ios_base::dec, std::ios_base::basefield);
    return io;
}

inline IOstream& hex(IOstream& io)
{
    io.setf(std::ios_base::hex, std::ios_base::basefield);
    return io;
}

inline IOstream& oct(IOstream& io)
{
    io.setf(std::ios_base::oct, std::ios_base::basefield);
    return io;
}

inline IOstream& fixed(IOstream& io)
{
    io.setf(std::ios_base::fixed, std::ios_base::floatfield);
    return io;
}

inline IOstream& scientific(IOstream& io)
{
    io.setf(std::ios_base::scientific, std::ios_base::floatfield);
    return io;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Detail
{

/*---------------------------------------------------------------------------*\
                   Class Detail::StreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A wrapper to hold a std::stream type for OpenFOAM wrapped streams.
//- This is necessary since the OpenFOAM streams hold a reference to
//- the normal std::stream
template<class StreamType>
class StreamAllocator
{
protected:

    // Protected Data

        //- The std::stream
        StreamType stream_;


    // Constructors

        //- Default construct (empty)
        StreamAllocator() = default;
};


} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Functions/Algorithms

namespace Detail
{

//- Termination for input looping (no-op)
template<class IS> inline void inputLoop(IS&) {}

//- Termination for output looping (no-op)
template<class OS> inline void outputLoop(OS&) {}

//- Input looping. Read into first parameter and recurse.
template<class IS, class Type, class... Args>
inline void inputLoop(IS& is, Type& arg1, Args&&... args)
{
    is >> arg1;
    Detail::inputLoop(is, std::forward<Args>(args)...);
}

//- Output looping. Write first parameter and recurse.
template<class OS, class Type, class... Args>
inline void outputLoop(OS& os, const Type& arg1, Args&&... args)
{
    os << arg1;
    Detail::outputLoop(os, std::forward<Args>(args)...);
}

} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
