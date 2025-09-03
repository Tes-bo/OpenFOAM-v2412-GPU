/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

InClass
    Foam::StringStream

Description
    Input/output from string buffers.

SourceFiles
    StringStream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_StringStream_H
#define Foam_StringStream_H

#include "ISstream.H"
#include "OSstream.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Input from string buffer, using a ISstream. Always UNCOMPRESSED.
class IStringStream
:
    public Foam::Detail::StreamAllocator<std::istringstream>,
    public Foam::ISstream
{
    typedef
        Foam::Detail::StreamAllocator<std::istringstream>
        allocator_type;

public:

    // Constructors

        //- Default construct or with specified stream option
        explicit IStringStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Construct from std::string
        explicit IStringStream
        (
            const std::string& s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {
            stream_.str(s);
        }

        //- Construct from char*
        explicit IStringStream
        (
            const char* s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {
            stream_.str(s);
        }

        //- Copy construct, copies content and format
        IStringStream(const IStringStream& str)
        :
            allocator_type(),
            ISstream(stream_, str.name(), static_cast<IOstreamOption>(str))
        {
            stream_.str(str.str());
        }


    // Member Functions

        //- Get the string.
        //- As Foam::string instead of std::string (may change in future)
        Foam::string str() const { return Foam::string(stream_.str()); }

        //- Set the string
        void str(const std::string& s) { stream_.str(s); }


        //- Reset the input buffer and rewind the stream
        virtual void reset(const std::string& s)
        {
            this->str(s);
            this->rewind();
        }

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const override;


    // Member Operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary.
        Istream& operator()() const
        {
            // Could also rewind
            return const_cast<IStringStream&>(*this);
        }


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct empty with given format
        explicit IStringStream(IOstreamOption::streamFormat fmt)
        :
            IStringStream(IOstreamOption(fmt))
        {}

        //- Construct from std::string with given format
        IStringStream
        (
            const std::string& s,
            IOstreamOption::streamFormat fmt
        )
        :
            IStringStream(s, IOstreamOption(fmt))
        {}

        //- Construct from char* with given format
        IStringStream
        (
            const char* s,
            IOstreamOption::streamFormat fmt
        )
        :
            IStringStream(s, IOstreamOption(fmt))
        {}

    #endif /* Foam_IOstream_extras */
};


/*---------------------------------------------------------------------------*\
                        Class OStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Output to string buffer, using a OSstream. Always UNCOMPRESSED.
class OStringStream
:
    public Foam::Detail::StreamAllocator<std::ostringstream>,
    public Foam::OSstream
{
    typedef
        Foam::Detail::StreamAllocator<std::ostringstream>
        allocator_type;

public:

    // Constructors

        //- Default construct or with specified stream option
        explicit OStringStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "output", streamOpt.format(), streamOpt.version())
        {}

        //- Copy construct, copies content and format
        OStringStream(const OStringStream& str)
        :
            allocator_type(),
            OSstream(stream_, str.name(), static_cast<IOstreamOption>(str))
        {
            stream_.str(str.str());
        }


    // Member Functions

        //- The number of bytes outputted
        std::streamsize count() { return stream_.tellp(); }

        //- Get the string.
        //- As Foam::string instead of std::string (may change in future)
        Foam::string str() const { return Foam::string(stream_.str()); }

        //- Set the string
        void str(const std::string& s) { stream_.str(s); }


        //- Reset the output buffer and rewind the stream
        void reset()
        {
            this->str("");  // No other way to reset the end
            this->rewind();
        }

        //- Rewind the output stream
        virtual void rewind()
        {
            stream_.rdbuf()->pubseekpos(0, std::ios_base::out);
        }

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const override;


    // Older style, without stream option (including 2012 release)
    #ifdef Foam_IOstream_extras

        //- Construct empty with given format
        explicit OStringStream(IOstreamOption::streamFormat fmt)
        :
            OStringStream(IOstreamOption(fmt))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
