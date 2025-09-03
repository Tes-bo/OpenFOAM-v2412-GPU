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
    Foam::ISpanStream

Description
    Similar to IStringStream but using an externally managed buffer for its
    input. This allows the input buffer to be filled (and refilled) from
    various sources.

    Note that this stream will normally be used as a "one-shot" reader.
    Caution must be exercised that the referenced buffer remains valid and
    without any intermediate resizing for the duration of the stream's use.

    An example of possible use:
    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        nread = something.read(buffer.data(),1024); // fill with content
        buffer.resize(nread);               // content size

        // Construct dictionary, or something else
        ISpanStream is(buffer)
        dictionary dict1(is);

        // Sometime later
        nread = something.read(buffer.data(),2048); // fill with content
        buffer.resize(nread);               // content size

        // Without intermediate variable
        dictionary dict2(ISpanStream(buffer)());
    \endcode

See Also
    Foam::ICharStream
    Foam::OCharStream
    Foam::OSpanStream

\*---------------------------------------------------------------------------*/

#ifndef Foam_ISpanStream_H
#define Foam_ISpanStream_H

#include "memoryStreamBuffer.H"
#include "ISstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ispanstream;
class ISpanStream;

// Older names (prior to 2023-08)
typedef ispanstream uiliststream;
typedef ISpanStream UIListStream;


/*---------------------------------------------------------------------------*\
                         Class ispanstream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::istringstream, but with an externally managed input buffer
//- which makes it most similar to std::ispanstream (C++23)
//  This allows the input buffer to be filled or refilled from various sources
//  without copying.
class ispanstream
:
    virtual public std::ios,
    protected Foam::memorybuf::in_base,
    public std::istream
{
    typedef Foam::memorybuf::in_base buffer_type;
    typedef std::istream stream_type;

public:

    // Constructors

        //- Default construct - empty
        ispanstream()
        :
            buffer_type(),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Construct (shallow copy) for character array and number of bytes
        ispanstream(const char* buffer, size_t nbytes)
        :
            buffer_type(const_cast<char*>(buffer), nbytes),
            stream_type(static_cast<buffer_type*>(this))
        {}

        #if __cplusplus >= 201703L
        //- Construct (shallow copy) from std::string_view content
        explicit ispanstream(std::string_view s)
        :
            buffer_type(const_cast<char*>(s.data()), s.size()),
            stream_type(static_cast<buffer_type*>(this))
        {}
        #endif

        //- Construct (shallow copy) from span character content
        explicit ispanstream(stdFoam::span<char> s)
        :
            buffer_type(const_cast<char*>(s.data()), s.size()),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Construct (shallow copy) from span character content
        explicit ispanstream(stdFoam::span<const char> s)
        :
            buffer_type(const_cast<char*>(s.data()), s.size()),
            stream_type(static_cast<buffer_type*>(this))
        {}


    // Member Functions

        //- The current get position within the buffer (tellg)
        std::streampos input_pos() const
        {
            return buffer_type::span_tellg();
        }

        //- The get buffer capacity
        std::streamsize capacity() const
        {
            return buffer_type::span_capacity();
        }

        //- The number of characters remaining in the get area.
        //- Same as (capacity() - input_pos())
        std::streamsize remaining() const
        {
            return buffer_type::span_remaining();
        }

        //- Span of the input characters (is modifiable!)
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

        //- For istringstream compatibility, return the buffer as string copy.
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
            buffer_type::pubseekpos(0, std::ios_base::in);
            stream_type::clear();  // Clear old errors
        }

        //- Reset the get buffer area
        void reset(const char* buffer, size_t nbytes)
        {
            buffer_type::resetg(const_cast<char*>(buffer), nbytes);
            stream_type::clear();  // Clear old errors
        }

        //- Reset the get buffer area to use the data from a string
        void reset(const std::string& s)
        {
            buffer_type::resetg(const_cast<char*>(&s[0]), s.size());
            stream_type::clear();  // Clear old errors
        }

        //- Some information about the input buffer position/capacity
        void debug_info(Ostream& os) const
        {
            os  << "get=" << input_pos() << '/' << capacity();
        }
};


/*---------------------------------------------------------------------------*\
                         Class ISpanStream Declaration
\*---------------------------------------------------------------------------*/

class ISpanStream
:
    public Foam::Detail::StreamAllocator<Foam::ispanstream>,
    public Foam::ISstream
{
    typedef
        Foam::Detail::StreamAllocator<Foam::ispanstream>
        allocator_type;

public:

    // Constructors

        //- Default construct (empty), optionally with specified stream option
        explicit ISpanStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Construct using specified buffer and number of bytes
        ISpanStream
        (
            const char* buffer,
            size_t nbytes,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ISpanStream(streamOpt)
        {
            reset(buffer, nbytes);
        }

        //- Use data area from string content
        explicit ISpanStream
        (
            const std::string& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ISpanStream(streamOpt)
        {
            reset(buffer);
        }

        //- Construct using data area from a List and its inherent storage size
        //  Uses addressed size, thus no special treatment for a DynamicList
        explicit ISpanStream
        (
            const ::Foam::UList<char>& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ISpanStream(buffer.cdata(), buffer.size(), streamOpt)
        {}

        #if __cplusplus >= 201703L
        //- Construct (shallow copy) from std::string_view content
        explicit ISpanStream
        (
            std::string_view s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ISpanStream(s.data(), s.size(), streamOpt)
        {}
        #endif

        //- Construct (shallow copy) from span character content
        explicit ISpanStream
        (
            stdFoam::span<const char> s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ISpanStream(s.data(), s.size(), streamOpt)
        {}

        //- Construct (shallow copy) from span character content
        explicit ISpanStream
        (
            stdFoam::span<char> s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ISpanStream(s.data(), s.size(), streamOpt)
        {}


    // Member Functions

        //- Position of the get buffer
        std::streampos tellg() const { return stream_.input_pos(); }

        //- The current get position within the buffer (tellg)
        std::streampos input_pos() const { return stream_.input_pos(); }

        //- The input list size. Same as capacity()
        label size() const { return label(stream_.capacity()); }

        //- The get buffer capacity
        std::streamsize capacity() const { return stream_.capacity(); }

        //- The number of characters remaining in the get area.
        //- Same as (capacity() - input_pos())
        std::streamsize remaining() const { return stream_.remaining(); }

        //- Span of the current input characters (is modifiable!)
        UList<char> list() const { return stream_.list(); }

        //- A string_view (c++17) or span view (older c++) of buffer contents
        auto view() const -> decltype(stream_.view())
        {
            return stream_.view();
        }

        //- For IStringStream compatibility, return the buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        auto str() const -> decltype(stream_.str())
        {
            return stream_.str();
        }

        //- Reset input area, position to buffer start and clear errors
        void reset(const char* buffer, size_t nbytes)
        {
            stream_.reset(buffer, nbytes);
            syncState();
        }

        //- Reset input area to use data from a string
        void reset(const std::string& s)
        {
            stream_.reset(s);
            syncState();
        }

        #if __cplusplus >= 201703L
        //- Reset input area to use data from a std::string_view
        void reset(std::string_view s)
        {
            stream_.reset(s.data(), s.size());
            syncState();
        }
        #endif

        //- Reset input area to use data from span character content
        void reset(stdFoam::span<char> s)
        {
            stream_.reset(s.data(), s.size());
            syncState();
        }

        //- Reset input area to use data from span character content
        void reset(stdFoam::span<const char> s)
        {
            stream_.reset(s.data(), s.size());
            syncState();
        }

        //- Rewind the stream, clearing any old errors
        virtual void rewind() override
        {
            stream_.rewind();
            syncState();
        }

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const override;


    // Member Operators

        //- A non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary
        Istream& operator()() const
        {
            // Could also rewind
            return const_cast<ISpanStream&>(*this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
