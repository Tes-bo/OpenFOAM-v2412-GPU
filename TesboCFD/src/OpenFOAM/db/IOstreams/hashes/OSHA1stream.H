/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::OSHA1stream

Description
    An output stream for calculating SHA1 digests.

\*---------------------------------------------------------------------------*/

#ifndef Foam_OSHA1stream_H
#define Foam_OSHA1stream_H

#include "SHA1.H"
#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class osha1stream Declaration
\*---------------------------------------------------------------------------*/

//- A basic output stream for calculating SHA1 digests
class osha1stream
:
    virtual public std::ios,
    public std::ostream
{
    //- A streambuf class for calculating SHA1 digests
    class sha1buf : public std::streambuf
    {
        //- This does all the work and has its own buffering
        SHA1 sha1_;

    protected:

        //- Output overflow handling - append to SHA1
        virtual int overflow(int_type c = traits_type::eof())
        {
            if (c != traits_type::eof()) sha1_.append(c);
            return c;
        }

        //- Put sequence of characters - append to SHA1
        virtual std::streamsize xsputn(const char* s, std::streamsize n)
        {
            if (n) sha1_.append(s, n);
            return n;
        }

    public:

        //- Default construct
        sha1buf() = default;

        //- Full access to the sha1
        SHA1& sha1() noexcept { return sha1_; }
    };


    // Private Data

        //- Reference to the underlying buffer
        sha1buf buf_;

public:

    // Constructors

        //- Default construct
        osha1stream() : std::ostream(&buf_) {}


    // Member Functions

        //- This hides both signatures of std::basic_ios::rdbuf()
        sha1buf* rdbuf() { return &buf_; }

        //- Full access to the sha1
        SHA1& sha1() noexcept { return buf_.sha1(); }

        //- Return SHA1::Digest for the data processed until now
        SHA1Digest digest() { return buf_.sha1().digest(); }

        //- Clear the SHA1 calculation
        void reset() { buf_.sha1().clear(); }
};


/*---------------------------------------------------------------------------*\
                         Class OSHA1stream Declaration
\*---------------------------------------------------------------------------*/

//- The output stream for calculating SHA1 digests
class OSHA1stream
:
    public Foam::Detail::StreamAllocator<Foam::osha1stream>,
    public Foam::OSstream
{
    typedef
        Foam::Detail::StreamAllocator<Foam::osha1stream>
        allocator_type;

public:

    //- No copy construct
    OSHA1stream(const OSHA1stream&) = delete;

    //- No copy assignment
    void operator=(const OSHA1stream&) = delete;


    // Constructors

        //- Construct with an empty digest
        explicit OSHA1stream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "sha1", streamOpt.format(), streamOpt.version())
        {}


    // Member Functions

        //- Full access to the sha1
        SHA1& sha1() noexcept { return stream_.sha1(); }

        //- Return SHA1::Digest for the data processed until now
        SHA1Digest digest() { return stream_.digest(); }

        //- Clear the SHA1 calculation
        void reset() { stream_.reset(); }


    // Write Functions

        //- Add (unquoted) string contents.
        //  Ensures that SHA1 of C-string or C++-string content are identical.
        virtual Ostream& write(const std::string& str) override
        {
            return writeQuoted(str, false);  // Unquoted!
        }


    // Housekeeping

        //- Deprecated(2017-07) clear the SHA1 calculation
        //  \deprecated(2017-07) - use reset() method
        void rewind()
        {
            stream_.sha1().clear();
        }


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct with an empty digest, using given format
        explicit OSHA1stream(IOstreamOption::streamFormat fmt)
        :
            OSHA1stream(IOstreamOption(fmt))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
