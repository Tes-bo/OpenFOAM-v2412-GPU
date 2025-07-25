/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::OTstream

Description
    A simple output token stream that can be used to build token lists.
    Always UNCOMPRESSED.

Note
    Appending single characters to token list is fragile.

SourceFiles
    OTstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OTstream_H
#define Foam_OTstream_H

#include "token.H"
#include "Ostream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OTstream Declaration
\*---------------------------------------------------------------------------*/

class OTstream
:
    public Ostream,
    public DynamicList<token>
{
public:

    // Constructors

        //- Default construct, set stream status
        explicit OTstream(IOstreamOption streamOpt = IOstreamOption())
        :
            Ostream(IOstreamOption(streamOpt.format(), streamOpt.version())),
            DynamicList<token>()
        {
            setOpened();
            setGood();
        }

        //- Copy construct
        OTstream(const OTstream& os)
        :
            Ostream(static_cast<IOstreamOption>(os)),
            DynamicList<token>(os.tokens())
        {
            setOpened();
            setGood();
        }

        //- Move construct
        OTstream(OTstream&& os)
        :
            Ostream(static_cast<IOstreamOption>(os)),
            DynamicList<token>(std::move(os.tokens()))
        {
            setOpened();
            setGood();
        }


    //- Destructor
    ~OTstream() = default;


    // Member Functions

        //- The tokens
        const DynamicList<token>& tokens() const noexcept { return *this; }

        //- The tokens
        DynamicList<token>& tokens() noexcept { return *this; }


    // Write

        //- Inherit write methods from Ostream
        using Ostream::writeQuoted;

        //- Write token to stream or otherwise handle it.
        //  \return false if the token type was not handled by this method
        virtual bool write(const token& tok) override;

        //- Write single character. Whitespace is suppressed.
        virtual Ostream& write(const char c) override;

        //- Write character/string content, with/without surrounding quotes
        virtual Ostream& writeQuoted
        (
            const char* str,
            std::streamsize len,
            const bool quoted=true
        ) override;

        //- Write the word-characters of a character string.
        //  Sends as a single char, or as word.
        virtual Ostream& write(const char* str) override;

        //- Write word
        virtual Ostream& write(const word& str) override;

        //- Write string
        virtual Ostream& write(const std::string& str) override;

        //- Write int32_t as a label
        virtual Ostream& write(const int32_t val) override;

        //- Write int64_t as a label
        virtual Ostream& write(const int64_t val) override;

        //- Write float
        virtual Ostream& write(const float val) override;

        //- Write double
        virtual Ostream& write(const double val) override;

        //- Write binary block with 8-byte alignment.
        virtual Ostream& write
        (
            const char* data,
            std::streamsize count
        ) override;

        //- Low-level raw binary output.
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
        virtual bool endRawWrite() override
        {
            return true;
        }

        //- Add indentation characters
        virtual void indent() override
        {}


    // Stream State Functions

        //- Return current stream flags.
        //- Dummy for token stream, returns 0.
        virtual std::ios_base::fmtflags flags() const override
        {
            return std::ios_base::fmtflags(0);
        }

        //- Set stream flags, return old stream flags.
        //- Dummy for token stream, returns 0.
        std::ios_base::fmtflags flags(std::ios_base::fmtflags) override
        {
            return std::ios_base::fmtflags(0);
        }

        //- Flush stream
        virtual void flush() override
        {}

        //- Add newline and flush stream
        virtual void endl() override
        {}

        //- Get the current padding character
        //  \return previous padding character
        virtual char fill() const override
        {
            return 0;
        }

        //- Set padding character for formatted field up to field width
        virtual char fill(const char) override
        {
            return 0;
        }

        //- Get width of output field
        virtual int width() const override
        {
            return 0;
        }

        //- Set width of output field
        //  \return previous width
        virtual int width(const int) override
        {
             return 0;
        }

        //- Get precision of output field
        virtual int precision() const override
        {
             return 0;
        }

        //- Set precision of output field
        //  \return old precision
        virtual int precision(const int) override
        {
             return 0;
        }


    // Other

        //- Rewind the output stream to position 0 (non-virtual!)
        //- and adjust the stream status (open/good/eof ...)
        //- Reset the output buffer and rewind the stream
        void reset() noexcept
        {
            DynamicList<token>::clear();
            setOpened();
            setGood();
        }

        //- Rewind the output stream to position 0
        //- and adjust the stream status (open/good/eof ...)
        virtual void rewind() { OTstream::reset(); }

        //- Print stream description to Ostream
        void print(Ostream& os) const override;


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct empty with format
        explicit OTstream(IOstreamOption::streamFormat fmt)
        :
            OTstream(IOstreamOption(fmt))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
