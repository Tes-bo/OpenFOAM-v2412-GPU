/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::OSpanStream

Description
    Similar to OStringStream but using an externally managed buffer for
    its output.

    This allows the output buffer to be reused and can make it easier when
    writing out data.  It is the user's responsibility to ensure proper
    synchronization in the sizes. Provided that the external buffer is large
    enough that overflow does not occur, the following usage pattern
    works.

    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        {
            OSpanStream os(buffer);
            os << "content1" << " and more content";
            buffer.resize(os.size());       // synchronize sizes
        }

        something.write(buffer, buffer.size());
    \endcode

    Although the OSpanStream is quite lightweight, there may be cases
    where it is preferable to reuse the stream as well.
    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        OSpanStream os(buffer);
        os << "content1" << " and more content";
        buffer.resize(os.size());           // synchronize sizes

        something.write(buffer, buffer.size());

        os.rewind();
        os << "content2";
        buffer.resize(os.size());      // synchronize sizes

        something.write(buffer, buffer.size());

        // or simply using the output size directly (without sync)
        os.rewind();
        os << "content3";

        something.write(buffer, os.size());
    \endcode

See Also
    Foam::ICharStream
    Foam::ISpanStream
    Foam::OCharStream

\*---------------------------------------------------------------------------*/

#ifndef Foam_OSpanStream_H
#define Foam_OSpanStream_H

#include "memoryStreamBuffer.H"
#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ospanstream;
class OSpanStream;

// Older names (prior to 2023-08)
typedef OSpanStream UOListStream;


/*---------------------------------------------------------------------------*\
                         Class ospanstream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::ostringstream, but with an externally managed output buffer
//- which makes it most similar to std::ospanstream (C++23)
class ospanstream
:
    virtual public std::ios,
    protected Foam::memorybuf::out_base,
    public std::ostream
{
    typedef Foam::memorybuf::out_base buffer_type;
    typedef std::ostream stream_type;

public:

    // Constructors

        //- Default construct - empty
        ospanstream()
        :
            buffer_type(),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Construct for character array and number of bytes
        ospanstream(char* buffer, size_t nbytes)
        :
            buffer_type(buffer, nbytes),
            stream_type(static_cast<buffer_type*>(this))
        {}


    // Member Functions

        //- The current output position within the buffer (tellp)
        std::streampos output_pos() const
        {
            return buffer_type::span_tellp();
        }

        //- The number of bytes outputted
        std::streamsize count() const
        {
            return buffer_type::size_bytes();
        }

        //- The put buffer capacity
        std::streamsize capacity() const
        {
            return buffer_type::span_capacity();
        }

        //- Span of the current output characters (is modifiable!)
        UList<char> list() const
        {
            return UList<char>
            (
                buffer_type::data_bytes(),
                label(buffer_type::size_bytes())
            );
        }

        //- A string_view (c++17) or span view (older c++) of buffer contents
        auto view() const -> decltype(buffer_type::view())
        {
            return buffer_type::view();
        }

        //- For ostringstream compatibility, return the buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        std::string str() const
        {
            return std::string
            (
                buffer_type::data_bytes(),
                buffer_type::size_bytes()
            );
        }

        //- Rewind the stream, clearing any old errors
        void rewind()
        {
            buffer_type::pubseekpos(0, std::ios_base::out);
            stream_type::clear();  // Clear any old errors
        }

        //- Reset the put buffer area
        void reset(char* buffer, size_t nbytes)
        {
            buffer_type::resetp(buffer, nbytes);
            stream_type::clear();  // Clear any old errors
        }

        //- Reset the put buffer area to use the data area from a string
        void reset(std::string& s)
        {
            s.resize(s.capacity());
            buffer_type::resetp(&s[0], s.size());
            stream_type::clear();  // Clear any old errors
        }

        //- Some information about the output buffer position/capacity
        void debug_info(Ostream& os) const
        {
            os  << "put=" << output_pos() << '/' << capacity();
        }
};


/*---------------------------------------------------------------------------*\
                         Class OSpanStream Declaration
\*---------------------------------------------------------------------------*/

class OSpanStream
:
    public Foam::Detail::StreamAllocator<Foam::ospanstream>,
    public Foam::OSstream
{
    typedef
        Foam::Detail::StreamAllocator<Foam::ospanstream>
        allocator_type;

public:

    // Constructors

        //- Default construct (empty output)
        explicit OSpanStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "output", streamOpt.format(), streamOpt.version())
        {}

        //- Use data area from string content
        explicit OSpanStream
        (
            std::string& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(streamOpt)
        {
            stream_.reset(buffer);
        }

        //- Construct using specified buffer and number of bytes
        OSpanStream
        (
            char* buffer,
            size_t nbytes,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(streamOpt)
        {
            stream_.reset(buffer, nbytes);
        }

        //- Construct using data area from a List and its inherent storage size
        explicit OSpanStream
        (
            ::Foam::UList<char>& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(buffer.data(), buffer.size(), streamOpt)
        {}

        //- Construct using full data area from DynamicList
        template<int SizeMin>
        explicit OSpanStream
        (
            ::Foam::DynamicList<char,SizeMin>& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(buffer.data(), buffer.capacity(), streamOpt)
        {
            buffer.resize(buffer.capacity());  // Uses entire space
        }


    // Member Functions

        //- Position of the put buffer
        std::streampos tellp() const { return stream_.output_pos(); }

        //- The current output position within the buffer (tellp)
        std::streampos output_pos() const { return stream_.output_pos(); }

        //- The number of bytes outputted
        std::streamsize count() const { return stream_.count(); }

        //- The current output size. Same as count(), output_pos(), tellp().
        label size() const { return label(stream_.count()); }

        //- The put buffer capacity
        std::streamsize capacity() const { return stream_.capacity(); }

        //- Span of the current output characters (is modifiable!)
        UList<char> list() const { return stream_.list(); }

        //- A string_view (c++17) or span view (older c++) of buffer contents
        auto view() const -> decltype(stream_.view())
        {
            return stream_.view();
        }

        //- For OStringStream compatibility, return buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        auto str() const -> decltype(stream_.str())
        {
            return stream_.str();
        }

        //- Reset the put area
        void reset(char* buffer, size_t nbytes)
        {
            stream_.reset(buffer, nbytes);
            syncState();
        }

        //- Reset the put buffer area to use the data area from a string
        void reset(std::string& s)
        {
            stream_.reset(s);
            syncState();
        }

        //- Rewind the stream, clearing any old errors
        virtual void rewind()
        {
            stream_.rewind();
            syncState();
        }

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
