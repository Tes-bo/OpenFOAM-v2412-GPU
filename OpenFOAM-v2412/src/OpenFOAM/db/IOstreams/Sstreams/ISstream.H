/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::ISstream

Description
    Generic input stream using a standard (STL) stream.

SourceFiles
    ISstreamI.H
    ISstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ISstream_H
#define Foam_ISstream_H

#include "Istream.H"
#include "fileName.H"
#include <limits>
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ISstream Declaration
\*---------------------------------------------------------------------------*/

class ISstream
:
    public Istream
{
    // Private Data

        //- The input stream path
        fileName name_;

        //- The input stream
        std::istream& is_;


    // Private Member Functions

        //- Get the next valid (non-whitespace) character,
        //- after skipping any C/C++ comments.
        char nextValid();

        //- Read into compound token (assumed to be a known type)
        virtual bool readCompoundToken(token& tok, const word& compoundType);

        //- No copy assignment
        void operator=(const ISstream&) = delete;


public:

    // Constructors

        //- Construct wrapper around std::istream, set stream status
        //  Default stream options (ASCII, uncompressed)
        inline ISstream
        (
            std::istream& is,
            const string& streamName,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct wrapper around std::istream, set stream status
        ISstream
        (
            std::istream& is,
            const string& streamName,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            ISstream(is, streamName, IOstreamOption(fmt, cmp))
        {}

        //- Construct wrapper around std::istream, set stream status
        ISstream
        (
            std::istream& is,
            const string& streamName,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            ISstream(is, streamName, IOstreamOption(fmt, ver, cmp))
        {}


    //- Destructor
    virtual ~ISstream() = default;


    // Member Functions

    // Characteristics

        //- The name of the input serial stream.
        //- (eg, the name of the Fstream file name)
        virtual const fileName& name() const override { return name_; }

        //- The name of the input serial stream, for modification.
        //  Use with caution since some classes (eg, Fstream)
        //  also use this for filesystem information!
        virtual fileName& name() { return name_; }


    // STL stream

        //- Const access to underlying std::istream
        virtual const std::istream& stdStream() const { return is_; }

        //- Access to underlying std::istream
        virtual std::istream& stdStream() { return is_; }


    // Stream State

        //- Return current stream flags
        virtual std::ios_base::fmtflags flags() const override
        {
            return is_.flags();
        }

        //- Set stream flags, return old stream flags
        virtual std::ios_base::fmtflags flags
        (
            std::ios_base::fmtflags f
        ) override
        {
            return is_.flags(f);
        }

        //- Set stream state to match that of the std::istream
        void syncState()
        {
            setState(is_.rdstate());
        }


    // Special-purpose Functions

        //- Discard until end of C-style comment '*/'
        //  \return False if stream exhausted before finding the comment end
        bool seekCommentEnd_Cstyle();

        //- Raw, low-level get into a string.
        //- Continues reading \b after an initial left-brace until it finds
        //- the matching closing right-brace.
        //  Tracks balanced pairs, trims out leading/trailing whitespace.
        //
        //  \return False if stream exhausted before finding closing brace
        bool continueReadUntilRightBrace
        (
            std::string& str,
            const bool stripComments = true
        );


    // Serial-stream functions

        //- Raw, low-level get character function.
        inline ISstream& get(char& c);

        //- Raw, low-level peek function.
        //  Does not remove the character from the stream.
        //  Returns the next character in the stream or EOF if the
        //  end of file is read.
        inline int peek();

        //- Raw, low-level getline (until delimiter) into a string.
        inline ISstream& getLine(std::string& str, char delim = '\n');

        //- Low-level discard until delimiter
        //  \return the number of characters extracted
        inline std::streamsize getLine(std::nullptr_t, char delim = '\n');

        //- Raw, low-level putback character function.
        inline ISstream& putback(const char c);


    // Read Functions

        //- Return next token from stream
        virtual Istream& read(token& t) override;

        //- Read a character
        virtual Istream& read(char& c) override;

        //- Read a word
        virtual Istream& read(word& str) override;

        //- Read a string (including enclosing double-quotes).
        //  Backslashes are retained, except when escaping double-quotes
        //  and an embedded newline character.
        virtual Istream& read(string& str) override;

        //- Read a label
        virtual Istream& read(label& val) override;

        //- Read a float
        virtual Istream& read(float& val) override;

        //- Read a double
        virtual Istream& read(double& val) override;

        //- Read binary block (with any possible block delimiters).
        //- Reading into a null pointer behaves like a forward seek of
        //- count characters.
        virtual Istream& read(char* data, std::streamsize count) override;

        //- Low-level raw binary read (without possible block delimiters).
        //- Reading into a null pointer behaves like a forward seek of
        //- count characters.
        virtual Istream& readRaw
        (
            char* data,
            std::streamsize count
        ) override;

        //- Start of low-level raw binary read
        virtual bool beginRawRead() override;

        //- End of low-level raw binary read
        virtual bool endRawRead() override;

        //- Rewind the stream so that it may be read again
        virtual void rewind() override;


    // Print

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
