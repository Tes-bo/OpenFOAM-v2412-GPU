/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::OSstream

Description
    Generic output stream using a standard (STL) stream.

SourceFiles
    OSstreamI.H
    OSstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OSstream_H
#define Foam_OSstream_H

#include "Ostream.H"
#include "fileName.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OSstream Declaration
\*---------------------------------------------------------------------------*/

class OSstream
:
    public Ostream
{
    // Private Data

        //- The output stream path
        fileName name_;

        //- The output stream
        std::ostream& os_;


public:

    // Generated Methods

        //- Copy construct
        OSstream(const OSstream&) = default;

        //- No copy assignment
        void operator=(const OSstream&) = delete;


    // Constructors

        //- Construct wrapper around std::ostream, set stream status
        //  Default stream options (ASCII, uncompressed)
        inline OSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct wrapper around std::ostream, set stream status
        OSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            OSstream(os, streamName, IOstreamOption(fmt, cmp))
        {}

        //- Construct wrapper around std::ostream, set stream status
        OSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            OSstream(os, streamName, IOstreamOption(fmt, ver, cmp))
        {}


    // Member Functions

    // Characteristics

        //- Get the name of the output serial stream.
        //- (eg, the name of the Fstream file name)
        virtual const fileName& name() const override { return name_; }

        //- The name of the output serial stream, for modification.
        //  Use with caution since some classes (eg, Fstream)
        //  also use this for filesystem information!
        virtual fileName& name() { return name_; }


    // STL stream

        //- Const access to underlying std::ostream
        virtual const std::ostream& stdStream() const { return os_; }

        //- Access to underlying std::ostream
        virtual std::ostream& stdStream() { return os_; }


    // Stream State

        //- Get current stream flags
        virtual std::ios_base::fmtflags flags() const override
        {
            return os_.flags();
        }

        //- Set stream flags, return old stream flags
        virtual std::ios_base::fmtflags flags
        (
            std::ios_base::fmtflags f
        ) override
        {
            return os_.flags(f);
        }

        //- Set stream state to match that of the std::ostream
        void syncState()
        {
            setState(os_.rdstate());
        }


    // Write Functions

        //- Inherit write methods from Ostream
        using Ostream::writeQuoted;

        //- Write token to stream or otherwise handle it.
        //  \return false if the token type was not handled by this method
        virtual bool write(const token& tok) override;

        //- Write character
        virtual Ostream& write(const char c) override;

        //- Write character/string content, with/without surrounding quotes
        virtual Ostream& writeQuoted
        (
            const char* str,
            std::streamsize len,
            const bool quoted=true
        ) override;

        //- Write character string
        virtual Ostream& write(const char* str) override;

        //- Write word
        virtual Ostream& write(const word& str) override;

        //- Write string (quoted)
        //  In the rare case that the string contains a final trailing
        //  backslash, it will be dropped to the appearance of an escaped
        //  double-quote.
        virtual Ostream& write(const std::string& str) override;

        //- Write int32_t
        virtual Ostream& write(const int32_t val) override;

        //- Write int64_t
        virtual Ostream& write(const int64_t val) override;

        //- Write float
        virtual Ostream& write(const float val) override;

        //- Write double
        virtual Ostream& write(const double val) override;

        //- Write binary block
        virtual Ostream& write
        (
            const char* data,
            std::streamsize count
        ) override;

        //- Low-level raw binary output
        virtual Ostream& writeRaw
        (
            const char* data,
            std::streamsize count
        ) override;

        //- Begin marker for low-level raw binary output.
        //  The count indicates the number of bytes for subsequent
        //  writeRaw calls.
        virtual bool beginRawWrite(std::streamsize count) override;

        //- End marker for low-level raw binary output.
        virtual bool endRawWrite() override;

        //- Add indentation characters
        virtual void indent() override;


    // Stream state functions

        //- Flush stream
        virtual void flush() override;

        //- Add newline and flush stream
        virtual void endl() override;

        //- Get the current padding character
        virtual char fill() const override;

        //- Set padding character for formatted field up to field width
        //  \return previous padding character
        virtual char fill(const char fillch) override;

        //- Get width of output field
        virtual int width() const override;

        //- Set width of output field
        //  \return previous width
        virtual int width(const int w) override;

        //- Get precision of output field
        virtual int precision() const override;

        //- Set precision of output field
        //  \return old precision
        virtual int precision(const int p) override;


    // Print

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OSstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
