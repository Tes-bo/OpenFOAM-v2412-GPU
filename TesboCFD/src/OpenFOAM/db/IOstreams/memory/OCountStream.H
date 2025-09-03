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
    Foam::OSCountStream

Description
    An output stream for calculating byte counts.

\*---------------------------------------------------------------------------*/

#ifndef Foam_OScountStream_H
#define Foam_OScountStream_H

#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ocountstream Declaration
\*---------------------------------------------------------------------------*/

//- Trivial output stream for calculating byte counts.
//  Since all output values are discarded, it can be used as a /dev/null
//  output buffer as well.
class ocountstream
:
    virtual public std::ios,
    public std::ostream
{
    //- A streambuf class for determining byte counts
    class countbuf : public std::streambuf
    {
        //- The number of bytes counted
        std::streamsize size_;

    protected:

        //- Set position pointer to relative position
        virtual std::streampos seekoff
        (
            std::streamoff off,
            std::ios_base::seekdir way,
            std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
        )
        {
            if (which & std::ios_base::out)
            {
                if (way == std::ios_base::beg)
                {
                    size_ = off;
                }
                else if (way == std::ios_base::cur)
                {
                    size_ += off;
                }
                else if (way == std::ios_base::end)
                {
                    // not really possible
                }

                return size_;  // Like span_tellp()
            }

            return -1;
        }

        //- Set position pointer to absolute position
        //  For the counter, adjust the count accordingly.
        virtual std::streampos seekpos
        (
            std::streampos pos,
            std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
        )
        {
            return seekoff(pos, std::ios_base::beg, which);
        }

        //- Output overflow handling - increment counter
        virtual int overflow(int_type c = traits_type::eof())
        {
            if (c != traits_type::eof()) ++size_;
            return c;
        }

        //- Put sequence of characters - increment counter
        virtual std::streamsize xsputn(const char* s, std::streamsize n)
        {
            size_ += n;
            return n;
        }

    public:

        // Constructors

            //- Default construct, count = 0
            countbuf() : size_(0) {}


        // Member Functions

            //- The number of bytes counted
            std::streamsize count() const noexcept { return size_; }

            //- Reset the count
            void reset(std::streamsize n = 0) noexcept { size_ = n; }
    };


    // Private Data

        typedef countbuf buffer_type;
        typedef std::ostream stream_type;

        //- Reference to the underlying buffer
        buffer_type buf_;

public:

    // Constructors

        //- Default construct
        ocountstream() : stream_type(&buf_) {}


    // Member Functions

        //- This hides both signatures of std::basic_ios::rdbuf()
        countbuf* rdbuf() { return &buf_; }

        //- The number of bytes counted
        std::streamsize count() const noexcept { return buf_.count(); }

        //- Reset the count
        void reset(std::streamsize n = 0) noexcept
        {
            buf_.reset(n);
            stream_type::clear();  // Clear old errors
        }

        //- Some information about the output buffer position/capacity
        void debug_info(Ostream& os) const
        {
            os  << "count=" << buf_.count();
        }
};


/*---------------------------------------------------------------------------*\
                        Class OCountStream Declaration
\*---------------------------------------------------------------------------*/

//- An output stream for calculating byte counts
class OCountStream
:
    public Foam::Detail::StreamAllocator<Foam::ocountstream>,
    public Foam::OSstream
{
    typedef
        Foam::Detail::StreamAllocator<Foam::ocountstream>
        allocator_type;

public:

    // Constructors

        //- Default construct
        explicit OCountStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "count", streamOpt.format(), streamOpt.version())
        {}

        //- Copy construct
        OCountStream(const OCountStream& str)
        :
            allocator_type(),
            OSstream(stream_, str.name(), static_cast<IOstreamOption>(str))
        {
            stream_.reset(str.count());
        }


    // Member Functions

        //- The number of bytes counted
        std::streamsize count() const noexcept { return stream_.count(); }

        //- The number of bytes counted
        std::streamsize size() const noexcept { return stream_.count(); }

        //- Reset the count
        void reset(std::streamsize n = 0) noexcept { stream_.reset(n); }

        //- Rewind the stream, reset the count, clearing any old errors
        virtual void rewind()
        {
            stream_.reset();
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
