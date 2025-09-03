/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::dummyISstream

Description
    Dummy input stream, which can be used as a placeholder for interfaces
    taking an Istream or ISstream. Aborts at any attempt to read from it.

Note
    The inheritance from an empty IStringStream is solely for convenience
    of implementation and should not be relied upon.

SourceFiles
    dummyISstream.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dummyISstream_H
#define Foam_dummyISstream_H

#include "StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dummyISstream Declaration
\*---------------------------------------------------------------------------*/

class dummyISstream
:
    public IStringStream
{
public:

    // Generated Methods

        //- Default construct
        dummyISstream() = default;

        //- Destructor
        virtual ~dummyISstream() = default;


    // Member Functions

    // Stream-state

        //- Return current stream flags.
        //- Dummy for dummy stream, returns 0.
        virtual std::ios_base::fmtflags flags() const override
        {
            return std::ios_base::fmtflags(0);
        }

        //- Set stream flags, return old stream flags.
        //- Dummy for dummy stream, returns 0.
        virtual std::ios_base::fmtflags flags(std::ios_base::fmtflags) override
        {
            return std::ios_base::fmtflags(0);
        }


    // Read Functions

        //- Return next token from stream
        virtual Istream& read(token&) override
        {
            NotImplemented;
            return *this;
        }

        //- Read a character
        virtual Istream& read(char&) override
        {
            NotImplemented;
            return *this;
        }

        //- Read a word
        virtual Istream& read(word&) override
        {
            NotImplemented;
            return *this;
        }

        // Read a string (including enclosing double-quotes)
        virtual Istream& read(string&) override
        {
            NotImplemented;
            return *this;
        }

        //- Read a label
        virtual Istream& read(label&) override
        {
            NotImplemented;
            return *this;
        }

        //- Read a float
        virtual Istream& read(float&) override
        {
            NotImplemented;
            return *this;
        }

        //- Read a double
        virtual Istream& read(double&) override
        {
            NotImplemented;
            return *this;
        }

        //- Read binary block
        virtual Istream& read(char*, std::streamsize) override
        {
            NotImplemented;
            return *this;
        }

        //- Low-level raw binary read
        virtual Istream& readRaw(char*, std::streamsize) override
        {
            NotImplemented;
            return *this;
        }

        //- Start of low-level raw binary read
        virtual bool beginRawRead() override
        {
            return false;
        }

        //- End of low-level raw binary read
        virtual bool endRawRead() override
        {
            return false;
        }

        //- Rewind the stream so that it may be read again
        virtual void rewind() override
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
