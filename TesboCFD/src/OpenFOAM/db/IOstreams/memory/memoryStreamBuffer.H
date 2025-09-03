/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::memorybuf

Description
    A std::streambuf used for memory buffer streams such as
    ispanstream, ocharstream, etc.

\*---------------------------------------------------------------------------*/

#ifndef Foam_memoryStreamBuffer_H
#define Foam_memoryStreamBuffer_H

#include "stdFoam.H"  // For span
#include "DynamicList.H"

#include <memory>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class memorybuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf for memory similar to std::spanbuf (C++23)
class memorybuf
:
    public std::streambuf
{
protected:

    //- Set position pointer to relative position
    virtual std::streampos seekoff
    (
        std::streamoff off,
        std::ios_base::seekdir way,
        std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
    )
    {
        const bool testin  = which & std::ios_base::in;
        const bool testout = which & std::ios_base::out;

        if (way == std::ios_base::beg)
        {
            if (testin)
            {
                setg(eback(), eback(), egptr());
                gbump(off);
            }
            if (testout)
            {
                setp(pbase(), epptr());
                pbump(off);
            }
        }
        else if (way == std::ios_base::cur)
        {
            if (testin)
            {
                gbump(off);
            }
            if (testout)
            {
                pbump(off);
            }
        }
        else if (way == std::ios_base::end)
        {
            if (testin)
            {
                setg(eback(), eback(), egptr());
                gbump(egptr() - eback() - off);
            }
            if (testout)
            {
                setp(pbase(), epptr());
                pbump(epptr() - pbase() - off);
            }
        }

        if (testin)
        {
            return (gptr() - eback()); // span_tellg()
        }
        if (testout)
        {
            return (pptr() - pbase()); // span_tellp()
        }

        return -1;
    }


    //- Set position pointer to absolute position
    virtual std::streampos seekpos
    (
        std::streampos pos,
        std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
    )
    {
        return seekoff(pos, std::ios_base::beg, which);
    }

public:

    // Forward Declarations
    class in_base;
    class in_dynamic;
    class out_base;
    class out_dynamic;
};


/*---------------------------------------------------------------------------*\
                     Class memorybuf::in_base Declaration
\*---------------------------------------------------------------------------*/

//- The base input streambuf with memory access
class memorybuf::in_base
:
    public memorybuf
{
protected:

    //- Get sequence of characters from a fixed region
    virtual std::streamsize xsgetn(char* s, std::streamsize n)
    {
        std::streamsize count = 0;
        while (count < n && gptr() < egptr())
        {
            *(s + count++) = *(gptr());
            gbump(1);
        }
        return count;
    }

public:

    // Constructors

        //- Default construct
        in_base() = default;

        //- Construct for character array (can be nullptr) and number of bytes
        in_base(char* s, std::streamsize n)
        {
            resetg(s, n);
        }


    // Member Functions

        // //- Reset get buffer pointer to the beginning of existing span
        // void rewind() { setg(eback(), eback(), egptr()); }

        //- Reset get buffer with character data (can be nullptr) and count
        //  Sets get pointer to the begin.
        void resetg(char* s, std::streamsize n)
        {
            if (s)
            {
                setg(s, s, s + n);
            }
            else
            {
                setg(nullptr, nullptr, nullptr);
            }
        }

        //- The current buffer get position
        std::streamsize span_tellg() const { return (gptr() - eback()); }

        //- The get buffer capacity
        std::streamsize span_capacity() const { return (egptr() - eback()); }

        //- The number of characters remaining in the get area
        std::streamsize span_remaining() const
        {
            return (gptr() < egptr()) ? (egptr() - gptr()) : 0;
        }

        //- The span data (start of input characters)
        char* data_bytes() const { return eback(); }

        //- The span size (number of input characters)
        std::streamsize size_bytes() const { return (egptr() - eback()); }

        #if __cplusplus >= 201703L
        std::string_view view() const
        {
            return std::string_view(data_bytes(), size_bytes());
        }
        #else
        stdFoam::span<const char> view() const
        {
            return stdFoam::span<const char>(data_bytes(), size_bytes());
        }
        #endif

        //- Some information about the input buffer position/capacity
        void info(Ostream& os) const
        {
            os  << "get=" << span_tellg() << '/' << span_capacity();
        }
};


/*---------------------------------------------------------------------------*\
                    Class memorybuf::in_dynamic Declaration
\*---------------------------------------------------------------------------*/

//- An output streambuf for memory access
class memorybuf::in_dynamic
:
    public memorybuf::in_base
{
private:

    //- Character storage
    List<char> storage_;


public:

    // Constructors

        //- Default construct - empty
        in_dynamic() = default;

        //- Copy construct from content
        in_dynamic(const char* s, std::streamsize n)
        {
            if (s && n)
            {
                storage_.resize_nocopy(n);
                std::copy(s, (s + n), storage_.data());
            }
            sync_gbuffer();
        }

        //- Move construct from List
        in_dynamic(::Foam::List<char>&& buffer)
        :
            storage_(std::move(buffer))
        {
            sync_gbuffer();
        }

        //- Move construct from DynamicList (added length only)
        template<int SizeMin>
        in_dynamic(::Foam::DynamicList<char,SizeMin>&& buffer)
        {
            storage_.transfer(buffer);  // Implies shrink_to_fit
            sync_gbuffer();
        }


    // Member Functions

        //- Sync get buffer pointers to agree with list dimensions
        //  Sets get pointer to the begin (rewind).
        void sync_gbuffer()
        {
            resetg(storage_.data(), storage_.size());
        }

        //- Reset content (copy)
        void reset(const char* s, std::streamsize n)
        {
            if (s && n)
            {
                storage_.resize_nocopy(n);
                std::copy(s, (s + n), storage_.data());
            }
            else
            {
                storage_.clear();
            }
            sync_gbuffer();
        }

        //- Exchange buffer content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            other.swap(storage_);  // Swap contents
            sync_gbuffer();
        }

        //- Exchange buffer content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            // NB: not storage_.swap(other)! - incorrect slicing
            other.swap(storage_);  // Swap contents: implies shrink_to_fit
            sync_gbuffer();
        }

        //- Reset buffer and return contents
        DynamicList<char> release()
        {
            DynamicList<char> chars(std::move(storage_));
            sync_gbuffer();
            return chars;
        }
};


/*---------------------------------------------------------------------------*\
                     Class memorybuf::out_base Declaration
\*---------------------------------------------------------------------------*/

//- An output streambuf for memory access
class memorybuf::out_base
:
    public memorybuf
{
protected:

    //- Put sequence of characters to a fixed region
    virtual std::streamsize xsputn(const char* s, std::streamsize n)
    {
        std::streamsize count = 0;
        while (count < n && pptr() < epptr())
        {
            *(pptr()) = *(s + count++);
            pbump(1);
        }
        return count;
    }

public:

    // Constructors

        //- Default construct
        out_base() = default;

        //- Construct for character array (can be nullptr) and number of bytes
        out_base(char* s, std::streamsize n)
        {
            resetp(s, n);
        }


    // Member Functions

        // //- Reset put buffer pointer to the beginning of existing span
        // void rewind() { setp(pbase(), epptr()); }

        //- Reset put buffer with character data (can be nullptr) and count
        //  Sets put pointer to the begin.
        inline void resetp(char* s, std::streamsize n)
        {
            if (s)
            {
                // As per (std::ios_base::out && !std::ios_base::ate)
                setp(s, s + n);
                // No treatment for (std::ios_base::out && std::ios_base::ate)
            }
            else
            {
                setp(nullptr, nullptr);
            }
        }

        //- The current buffer put position
        std::streamsize span_tellp() const { return (pptr() - pbase()); }

        //- The put buffer capacity
        std::streamsize span_capacity() const { return (epptr() - pbase()); }

        //- The span data (start of output characters)
        char* data_bytes() const { return pbase(); }

        //- The span size (size of output buffer)
        std::streamsize size_bytes() const { return (pptr() - pbase()); }

        #if __cplusplus >= 201703L
        std::string_view view() const
        {
            return std::string_view(data_bytes(), size_bytes());
        }
        #else
        stdFoam::span<const char> view() const
        {
            return stdFoam::span<const char>(data_bytes(), size_bytes());
        }
        #endif

        //- Some information about the output buffer position/capacity
        void info(Ostream& os) const
        {
            os  << "put=" << span_tellp() << '/' << span_capacity();
        }
};


/*---------------------------------------------------------------------------*\
                   Class memorybuf::out_dynamic Declaration
\*---------------------------------------------------------------------------*/

//- An output streambuf for memory access
class memorybuf::out_dynamic
:
    public memorybuf::out_base
{
private:

    //- Character storage.
    //  Internally manage like a DynamicList, with its capacity known
    //  from the list size and the addressable size known through the
    //  stream pointers.
    List<char> storage_;

protected:

    //- Handle overflow
    virtual int overflow(int_type c = traits_type::eof())
    {
        if (c != traits_type::eof())
        {
            // Need more space?
            reserve(1 + span_tellp());

            *(pptr()) = c;
            pbump(1);
        }
        return c;
    }

    //- Put sequence of characters
    virtual std::streamsize xsputn(const char* s, std::streamsize n)
    {
        // Enough space so that appends work without problem
        reserve(n + span_tellp());

        std::streamsize count = 0;
        while (count < n && pptr() < epptr())
        {
            *(pptr()) = *(s + count++);
            pbump(1);
        }

        return count;
    }


public:

    // Constructors

        //- Default construct with initial reserved number of bytes.
        //  The value of 512 is a bit arbitrary, but consistent with
        //  std::stringstream
        out_dynamic(size_t nbytes = 512)
        :
            storage_(label(nbytes))
        {
            sync_pbuffer();
        }

        //- Move construct from List
        out_dynamic(::Foam::List<char>&& buffer)
        :
            storage_(std::move(buffer))
        {
            sync_pbuffer();
        }

        //- Move construct from DynamicList (uses entire capacity)
        template<int SizeMin>
        out_dynamic(::Foam::DynamicList<char,SizeMin>&& buffer)
        {
            buffer.resize(buffer.capacity());  // Use entire space
            storage_.transfer(buffer);
            sync_pbuffer();
        }


    // Member Functions

        //- Increment capacity (if needed) and adjust buffer pointers
        void reserve(const std::streamsize len)
        {
            if (storage_.size() < len)
            {
                const auto cur = span_tellp();   // Current location

                label newCapacity = 512;

                if (newCapacity < len)
                {
                    // Increase capacity (doubling)
                    newCapacity =
                        Foam::max(label(len), label(2*storage_.size()));
                }

                // Info<<"request:" << len
                //     << " cur cap:" << storage_.size()
                //     << " new cap:" << newCapacity
                //     << " pos:" << cur << endl;

                storage_.resize(newCapacity);
                sync_pbuffer();
                pbump(cur);
            }
        }

        //- Sync put buffer pointers to agree with list dimensions
        //  Sets put pointer to the begin (rewind).
        void sync_pbuffer()
        {
            resetp(storage_.data(), storage_.size());
        }

        //- Clear storage
        void clearStorage()
        {
            storage_.clear();
            sync_pbuffer();
        }

        //- Shrink storage to addressed storage
        void shrink()
        {
            const auto cur = span_tellp();  // Addressed length

            storage_.resize(cur);
            sync_pbuffer();
            pbump(cur);
        }

        //- Exchange buffer content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            const auto cur = span_tellp();  // Output length
            other.swap(storage_);
            other.resize(cur);              // Truncate to output length
            sync_pbuffer();
        }

        //- Exchange buffer content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            const auto cur = span_tellp();  // Output length

            other.resize(other.capacity()); // Use entire space
            other.swap(storage_);  // NB: not storage_.swap(other)
            other.resize(cur);              // Restrict to output length
            sync_pbuffer();
        }

        //- Reset buffer and return contents as a DynamicList.
        //- The list size corresponds to the region of output.
        DynamicList<char> release()
        {
            const auto cur = span_tellp();  // Output length
            DynamicList<char> chars(std::move(storage_));
            chars.resize(cur);              // Restrict to output length

            if (chars.empty()) chars.clearStorage();  // Can destroy now
            sync_pbuffer();
            return chars;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
