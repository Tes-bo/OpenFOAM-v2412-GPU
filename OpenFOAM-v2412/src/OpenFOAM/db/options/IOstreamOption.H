/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::IOstreamOption

Description
    A simple container for options an IOstream can normally have.

    The format (ASCII | BINARY) is typically controlled by enumerated
    names (ascii, binary).

    The compression (UNCOMPRESSED | COMPRESSED) is typically controlled
    by switch values (true/false, on/off, ...).

SourceFiles
    IOstreamOption.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOstreamOption_H
#define Foam_IOstreamOption_H

#include "word.H"
#include <ios>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class token;
class dictionary;
template<class EnumType> class Enum;

/*---------------------------------------------------------------------------*\
                       Class IOstreamOption Declaration
\*---------------------------------------------------------------------------*/

class IOstreamOption
{
public:

    // Public Data Types

        //- Data format (ascii | binary)
        enum streamFormat : char
        {
            ASCII = 0,          //!< "ascii" (normal default)
            BINARY              //!< "binary"
        };

        //- Compression treatment (UNCOMPRESSED | COMPRESSED)
        enum compressionType : char
        {
            UNCOMPRESSED = 0,   //!< compression = false
            COMPRESSED          //!< compression = true
        };

        //- File appending (NO_APPEND | APPEND_APP | APPEND_ATE)
        enum appendType : char
        {
            NO_APPEND = 0,      //!< no append (truncates existing)
            APPEND_APP,         //!< append (seek end each write)
            APPEND_ATE,         //!< append (seek end after open)
            NON_APPEND = NO_APPEND,     //!< old name for NO_APPEND
            APPEND = APPEND_APP         //!< old name for APPEND_APP
        };

        //- Atomic operations (output)
        enum atomicType : char
        {
            NON_ATOMIC = 0,     //!< atomic = false
            ATOMIC              //!< atomic = true
        };

        //- Float formats (eg, time directory name formats)
        enum class floatFormat : unsigned
        {
            //! default float notation
            general    = unsigned(0),
            //! fixed-point notation
            fixed      = unsigned(std::ios_base::fixed),
            //! scientific notation
            scientific = unsigned(std::ios_base::scientific)
        };


        //- Representation of a major/minor version number
        class versionNumber
        {
            //- The combined major/version number.
            short number_;

        public:

        // Constructors

            //- Default construct \em current version.
            //- The value (2.0) corresponds to the \em current version.
            constexpr versionNumber() noexcept
            :
                number_(20)
            {}

            //- Construct from major, number
            constexpr versionNumber(int major, int minor) noexcept
            :
                number_(10*major + (minor % 10))
            {}

            //- Construct from floating-point version number
            explicit constexpr versionNumber(const float ver) noexcept
            :
                number_(10*ver + 0.001) // Allow some rounding
            {}

            //- Construct by parsing string "major.minor"
            explicit versionNumber(const std::string& verNum);

            //- Construct from token (float, word, string)
            explicit versionNumber(const token& tok);

            //- Failsafe construct from dictionary lookup.
            versionNumber(const word& key, const dictionary& dict);


        // Member Functions

            //- A string representation as major.minor
            std::string str() const
            {
                return
                (
                    std::to_string(int(number_ / 10)) // major
                  + '.'
                  + std::to_string(int(number_ % 10)) // minor
                );
            }

            //- From version to canonical integer value
            int canonical() const noexcept
            {
                return number_;
            }

            //- From canonical integer value to version
            static versionNumber canonical(int verNum) noexcept
            {
                // Split into major/minor
                return versionNumber(int(verNum / 10), int(verNum % 10));
            }

            //- Compare differences in the versions
            //  Negative when 'this' is less than other.
            //  Positive when 'this' is greater than other.
            int compare(const versionNumber& other) const noexcept
            {
                return (number_ - other.number_);
            }
        };


    // Public Static Data

        //- Names for float formats (general, fixed, scientific)
        static const Enum<floatFormat> floatFormatNames;

        //- Stream format names (ascii, binary)
        static const Enum<streamFormat> formatNames;

        //- The current version number (2.0)
        static const versionNumber currentVersion;


    // Static Helpers

        //- Lookup floatFormat enum corresponding to the string
        //- (general | fixed | scientific).
        //
        //  If the string is not recognized, emit warning and return default.
        //  Silent if the string itself is empty.
        //
        //  \note Can be used as constructor substitute for the enumeration
        static floatFormat floatFormatEnum
        (
            const word& fmtName,
            const floatFormat deflt = floatFormat::general
        );

        //- getOrDefault floatFormat from dictionary,
        //- warn only on bad enumeration.
        static floatFormat floatFormatEnum
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const floatFormat deflt = floatFormat::general
        );

        //- Lookup streamFormat enum corresponding to the string
        //- (ascii | binary).
        //
        //  If the string is not recognized, emit warning and return default.
        //  Silent if the string itself is empty.
        //
        //  \note Can be used as constructor substitute for the enumeration
        static streamFormat formatEnum
        (
            const word& fmtName,
            const streamFormat deflt = streamFormat::ASCII
        );

        //- getOrDefault streamFormat from dictionary,
        //- warn only on bad enumeration.
        static streamFormat formatEnum
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const streamFormat deflt = streamFormat::ASCII
        );

        //- The compression enum corresponding to the string.
        //  Expects switch values (true/false, on/off, ...)
        //
        //  If the string is not recognized, emit warning and return default.
        //  Silent if the string itself is empty.
        //
        //  \note Can be used as constructor substitute for the enumeration
        static compressionType compressionEnum
        (
            const word& compName,
            const compressionType deflt = compressionType::UNCOMPRESSED
        );

        //- getOrDefault compressionType from dictionary,
        //- warn only on bad enumeration.
        static compressionType compressionEnum
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const compressionType deflt = compressionType::UNCOMPRESSED
        );


private:

    // Private Data

    // NB: ordered with versionNumber first (short) and
    // adjacent enums to minimize gaps

        //- Stream version number (eg, 2.0 for current dictionary format)
        versionNumber version_;

        //- Format: (ascii | binary)
        streamFormat format_;

        //- Compression: (on | off)
        compressionType compression_;


public:

    // Constructors

        //- Default construct (ASCII, UNCOMPRESSED, currentVersion)
        //- or construct with format, compression
        //  \note non-explicit for convenient construction
        constexpr IOstreamOption
        (
            streamFormat fmt = streamFormat::ASCII,
            compressionType comp = compressionType::UNCOMPRESSED
        ) noexcept
        :
            version_(),
            format_(fmt),
            compression_(comp)
        {}

        //- Construct from components (format, compression, version)
        constexpr IOstreamOption
        (
            streamFormat fmt,
            compressionType comp,
            versionNumber ver
        ) noexcept
        :
            version_(ver),
            format_(fmt),
            compression_(comp)
        {}

        //- Construct from components (format, version, compression)
        constexpr IOstreamOption
        (
            streamFormat fmt,
            versionNumber ver,
            compressionType comp = compressionType::UNCOMPRESSED
        ) noexcept
        :
            version_(ver),
            format_(fmt),
            compression_(comp)
        {}

        //- Copy construct with change of format
        IOstreamOption(const IOstreamOption& opt, streamFormat fmt) noexcept
        :
            version_(opt.version_),
            format_(fmt),
            compression_(opt.compression_)
        {}


    // Member Functions

        //- Get the current stream format
        streamFormat format() const noexcept { return format_; }

        //- Set the stream format
        //  \return the previous value
        streamFormat format(const streamFormat fmt) noexcept
        {
            streamFormat old(format_);
            format_ = fmt;
            return old;
        }

        //- Set the stream format from string value.
        //  If the string is not recognized, emit warning and leave unchanged.
        //  Silent if the string itself is empty.
        //  \return the previous value
        streamFormat format(const word& formatName)
        {
            streamFormat old(format_);
            format_ = formatEnum(formatName, format_);
            return old;
        }

        //- Get the stream compression
        compressionType compression() const noexcept { return compression_; }

        //- Set the stream compression
        //  \return the previous value
        compressionType compression(const compressionType comp) noexcept
        {
            compressionType old(compression_);
            compression_ = comp;
            return old;
        }

        //- Set the stream compression from string value.
        //  If the string is not recognized, emit warning and leave unchanged.
        //  Silent if the string itself is empty.
        //  \return the previous value
        compressionType compression(const word& compName)
        {
            compressionType old(compression_);
            compression_ = compressionEnum(compName, compression_);
            return old;
        }

        //- Get the stream version
        versionNumber version() const noexcept { return version_; }

        //- Set the stream version
        //  \return the previous value
        versionNumber version(const versionNumber ver) noexcept
        {
            versionNumber old(version_);
            version_ = ver;
            return old;
        }

        //- Set the stream version from token
        //  \return the previous value
        versionNumber version(const token& tok)
        {
            versionNumber old(version_);
            version_ = versionNumber(tok);
            return old;
        }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Output format type as text string (ascii | binary)
Ostream& operator<<(Ostream& os, const IOstreamOption::streamFormat& fmt);

//- Output version as major.minor text string
Ostream& operator<<(Ostream& os, const IOstreamOption::versionNumber& ver);


// * * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * //

// Comparison Operators

//- Version number equality
inline bool operator==
(
    const IOstreamOption::versionNumber& a,
    const IOstreamOption::versionNumber& b
) noexcept
{
    return a.compare(b) == 0;
}

//- Version number inequality
inline bool operator!=
(
    const IOstreamOption::versionNumber& a,
    const IOstreamOption::versionNumber& b
) noexcept
{
    return a.compare(b) != 0;
}

//- Version A older than B
inline bool operator<
(
    const IOstreamOption::versionNumber& a,
    const IOstreamOption::versionNumber& b
) noexcept
{
    return a.compare(b) < 0;
}

//- Version A same or older than B
inline bool operator<=
(
    const IOstreamOption::versionNumber& a,
    const IOstreamOption::versionNumber& b
) noexcept
{
    return a.compare(b) <= 0;
}

//- Version A newer than B
inline bool operator>
(
    const IOstreamOption::versionNumber& a,
    const IOstreamOption::versionNumber& b
) noexcept
{
    return a.compare(b) > 0;
}

//- Version A same or newer than B
inline bool operator>=
(
    const IOstreamOption::versionNumber& a,
    const IOstreamOption::versionNumber& b
) noexcept
{
    return a.compare(b) >= 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
