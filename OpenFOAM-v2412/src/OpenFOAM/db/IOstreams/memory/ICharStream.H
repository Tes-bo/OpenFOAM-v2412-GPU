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

Class
    Foam::ICharStream

Description
    An input stream that reads from a List and manages the List storage.
    Similar to IStringStream but with a List for its storage instead of
    as string to allow reuse of List contents without copying.

See Also
    Foam::OCharStream
    Foam::ISpanStream
    Foam::OSpanStream

\*---------------------------------------------------------------------------*/

#ifndef Foam_ICharStream_H
#define Foam_ICharStream_H

#include "ISpanStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class icharstream;
class ICharStream;
class OCharStream;

// Older names (prior to 2023-08)
typedef ICharStream IListStream;


/*---------------------------------------------------------------------------*\
                        Class icharstream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::istringstream, but with the ability to swap
//- character content.
//- Has some similarity to std::ispanstream (C++23)
class icharstream
:
    virtual public std::ios,
    protected Foam::memorybuf::in_dynamic,
    public std::istream
{
    typedef Foam::memorybuf::in_dynamic buffer_type;
    typedef std::istream stream_type;

public:

    // Constructors

        //- Default construct - empty
        icharstream()
        :
            buffer_type(),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Copy construct from content
        icharstream(const char* buffer, size_t nbytes)
        :
            icharstream()
        {
            reset(buffer, nbytes);
        }

        //- Move construct from List
        icharstream(List<char>&& buffer)
        :
            icharstream()
        {
            swap(buffer);
        }

        //- Move construct from DynamicList
        template<int SizeMin>
        icharstream(DynamicList<char,SizeMin>&& buffer)
        :
            icharstream()
        {
            swap(buffer);
        }


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

        //- Reset stream content (copy), reset positions
        void reset(const char* buffer, size_t nbytes)
        {
            buffer_type::reset(buffer, nbytes);
            stream_type::clear();  // Clear old errors
        }

        //- Exchange stream content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            buffer_type::swap(other);
            stream_type::clear();  // Clear old errors
        }

        //- Exchange stream content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            buffer_type::swap(other);
            stream_type::clear();  // Clear old errors
        }

        //- Reset stream and return contents as a List
        DynamicList<char> release()
        {
            DynamicList<char> chars(buffer_type::release());
            stream_type::clear();  // Clear old errors
            return chars;
        }

        //- Some information about the input buffer position/capacity
        void debug_info(Ostream& os) const
        {
            os  << "get=" << input_pos() << '/' << capacity();
        }
};


/*---------------------------------------------------------------------------*\
                         Class ICharStream Declaration
\*---------------------------------------------------------------------------*/

//- An ISstream with internal List storage. Always UNCOMPRESSED.
class ICharStream
:
    public Foam::Detail::StreamAllocator<Foam::icharstream>,
    public Foam::ISstream
{
    typedef
        Foam::Detail::StreamAllocator<Foam::icharstream>
        allocator_type;

public:

    // Constructors

        //- Default construct (empty), optionally with specified stream option
        explicit ICharStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Copy construct from string content
        explicit ICharStream
        (
            const std::string& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ICharStream(streamOpt)
        {
            stream_.reset(buffer.data(), buffer.size());
        }

        //- Move construct from List
        explicit ICharStream
        (
            List<char>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ICharStream(streamOpt)
        {
            stream_.swap(buffer);
        }

        //- Move construct from DynamicList (uses current size)
        template<int SizeMin>
        explicit ICharStream
        (
            DynamicList<char,SizeMin>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ICharStream(streamOpt)
        {
            stream_.swap(buffer);
        }


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

        //- Span of the input characters (is modifiable!)
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

        //- Reset content (copy)
        void reset(const char* buffer, size_t nbytes)
        {
            stream_.reset(buffer, nbytes);
            syncState();
        }

        //- Exchange stream content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            stream_.swap(other);
            syncState();
        }

        //- Exchange stream content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            stream_.swap(other);
            syncState();
        }

        //- Reset stream and return contents as a List
        DynamicList<char> release()
        {
            DynamicList<char> chars(stream_.release());
            syncState();
            return chars;
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
            return const_cast<ICharStream&>(*this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
