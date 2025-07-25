/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::Istream

Description
    An Istream is an abstract base class for all input systems
    (streams, files, token lists etc).  The basic operations
    are construct, close, read token, read primitive and read binary
    block.

    In addition, version control and line number counting is incorporated.
    Usually one would use the read primitive member functions, but if one
    were reading a stream on unknown data sequence one can read token by
    token, and then analyse.

SourceFiles
    Istream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Istream_H
#define Foam_Istream_H

#include "IOstream.H"
#include "token.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Istream Declaration
\*---------------------------------------------------------------------------*/

class Istream
:
    public IOstream
{
    // Private Data

        //- The last token put back on the stream
        token putBackToken_;

        //- Is a put-back token available?
        bool putBackAvail_;


protected:

    // Protected Member Functions

        //- True if putback token is in use
        bool hasPutback() const noexcept { return putBackAvail_; }


public:

    // Generated Methods

        //- Copy construct
        Istream(const Istream&) = default;

        //- Destructor
        virtual ~Istream() = default;


    // Constructors

        //- Default construct (ASCII, uncompressed),
        //- construct with specified stream option.
        explicit Istream(IOstreamOption streamOpt = IOstreamOption())
        :
            IOstream(streamOpt),
            putBackAvail_(false)
        {}

        //- Construct with format (uncompressed)
        explicit Istream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            Istream(IOstreamOption(fmt, cmp))
        {}


    // Member Functions

        // Token put-back

            //- Examine putback token without removing it.
            //  Returns const reference to \c token::undefinedToken
            //  if a putback is unavailable.
            const token& peekBack() const noexcept;

            //- Drop the putback token
            void putBackClear();

            //- Put back a token (copy). Only a single put back is permitted
            void putBack(const token& tok);

            //- Put back a token (move). Only a single put back is permitted
            void putBack(token&& tok);

            //- Retrieve the put-back token if there is one.
            //  \return false and sets token to undefined if no put-back
            //  was available
            bool getBack(token& tok);


       // Read Functions

            //- Return next token from stream
            virtual Istream& read(token&) = 0;

            //- Read a character
            virtual Istream& read(char&) = 0;

            //- Read a word
            virtual Istream& read(word&) = 0;

            //- Read a string (including enclosing double-quotes)
            virtual Istream& read(string&) = 0;

            //- Read a label
            virtual Istream& read(label&) = 0;

            //- Read a float
            virtual Istream& read(float&) = 0;

            //- Read a double
            virtual Istream& read(double&) = 0;

            //- Read binary block (with any possible block delimiters).
            //- Reading into a null pointer shall ideally behave like a seek
            //- operation.
            virtual Istream& read(char* data, std::streamsize count) = 0;

            //- Start of low-level raw binary read
            virtual bool beginRawRead() = 0;

            //- End of low-level raw binary read
            virtual bool endRawRead() = 0;

            //- Low-level raw binary read (without possible block delimiters).
            //- Reading into a null pointer shall ideally behave like a seek
            //- operation.
            virtual Istream& readRaw(char* data, std::streamsize count) = 0;

            //- Rewind the stream so that it may be read again
            virtual void rewind() = 0;


        // Read List punctuation tokens

            //- Begin read of data chunk, starts with '('.
            //  \return true or FatalIOError
            bool readBegin(const char* funcName);

            //- End read of data chunk, ends with ')'
            //  \return true or FatalIOError
            bool readEnd(const char* funcName);

            //- Begin read of list data, starts with '(' or '{'
            //  \return starting delimiter or FatalIOError
            char readBeginList(const char* funcName);

            //- End read of list data, ends with ')' or '}'
            //  \return closing delimiter or FatalIOError
            char readEndList(const char* funcName);


    // Member Operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary
        Istream& operator()() const;
};


// --------------------------------------------------------------------
// ------ Manipulators (not taking arguments)
// --------------------------------------------------------------------

//- An Istream manipulator
typedef Istream& (*IstreamManip)(Istream&);

//- operator>> handling for manipulators without arguments
inline Istream& operator>>(Istream& is, IstreamManip f)
{
    return f(is);
}

//- operator>> handling for manipulators without arguments
inline Istream& operator>>(Istream& is, IOstreamManip f)
{
    f(is);
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Detail
{
    //- Read binary block of contiguous data, possibly with conversion
    template<class T>
    void readContiguous(Istream& is, char* data, std::streamsize byteCount)
    {
        is.beginRawRead();

        if (is_contiguous_label<T>::value)
        {
            readRawLabel
            (
                is,
                reinterpret_cast<label*>(data),
                byteCount/sizeof(label)
            );
        }
        else if (is_contiguous_scalar<T>::value)
        {
            readRawScalar
            (
                is,
                reinterpret_cast<scalar*>(data),
                byteCount/sizeof(scalar)
            );
        }
        else
        {
            is.readRaw(data, byteCount);
        }

        is.endRawRead();
    }

} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Previously excluded (via token.H)
#ifdef NoRepository
    #include "HashTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
