/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::SHA1

Description
    Functions to compute SHA1 message digest according to the NIST
    specification FIPS-180-1.

    Adapted from the gnulib implementation.

See also
    Foam::SHA1Digest

SourceFiles
    SHA1I.H
    SHA1.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SHA1_H
#define Foam_SHA1_H

#include <string>
#include <cstdint>
#include "SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class SHA1 Declaration
\*---------------------------------------------------------------------------*/

class SHA1
{
    // Private Data

        //- Track if the hashsum has been finalized (added count, etc)
        bool finalized_;

        //- The hash sums
        uint32_t hashsumA_;
        uint32_t hashsumB_;
        uint32_t hashsumC_;
        uint32_t hashsumD_;
        uint32_t hashsumE_;

        //- The total number processed, saved as 64-bit
        uint32_t bufTotal_[2];

        //- The number of elements pending in the buffer
        uint32_t bufLen_;

        //- The input processing buffer
        uint32_t buffer_[32];


    // Private Member Functions

        //- Process data block-wise, LEN must be a multiple of 64!
        void processBlock(const void *data, size_t len);

        //- Process for the next LEN bytes, LEN need not be a multiple of 64.
        void processBytes(const void *data, size_t len);

        //- Calculate digest from current data.
        void calcDigest(SHA1Digest& dig) const;


public:

    // Constructors

        //- Default construct
        inline SHA1();

        //- Default construct and append initial string content
        inline explicit SHA1(const char* str);

        //- Default construct and append initial string content
        inline explicit SHA1(const std::string& str);


    // Member Functions

        //- Reset the hashed data before appending more
        void clear() noexcept;

        //- Append single character
        inline void append(char c);

        //- Append string for processing
        inline void append(const char* str);

        //- Append data for processing
        inline void append(const char* data, size_t len);

        //- Append string for processing
        inline void append(const std::string& str);

        //- Append substring for processing
        inline void append
        (
            const std::string& str,
            size_t pos,
            size_t len = std::string::npos
        );

        //- Finalized the calculations (normally not needed directly).
        //  Returns false if no bytes were passed for processing
        bool finalize();

        //- Calculate digest from current data.
        inline SHA1Digest digest() const;

        //- The digest (40-byte) text representation, optionally with '_' prefix
        inline std::string str(const bool prefixed=false) const;

        //- Write digest (40-byte) representation, optionally with '_' prefix
        inline Ostream& write(Ostream& os, const bool prefixed=false) const;


    // Member Operators

        //- Cast conversion to a SHA1Digest,
        //  calculates the %digest from the current data
        inline operator SHA1Digest() const;

        //- Equality operator, compares %digests
        inline bool operator==(const SHA1& rhs) const;

        //- Compare %digest
        inline bool operator==(const SHA1Digest& dig) const;

        //- Compare %digest to (40-byte) text representation (eg, from sha1sum)
        //  An %empty string is equivalent to
        //  "0000000000000000000000000000000000000000"
        inline bool operator==(const std::string& hexdigits) const;

        //- Compare %digest to (40-byte) text representation (eg, from sha1sum)
        //  A %null or %empty string is equivalent to
        //  "0000000000000000000000000000000000000000"
        inline bool operator==(const char* hexdigits) const;

        //- Inequality operator, compares %digests
        inline bool operator!=(const SHA1&) const;

        //- Inequality operator, compare %digest
        inline bool operator!=(const SHA1Digest&) const;

        //- Inequality operator, compares %digests
        inline bool operator!=(const std::string& hexdigits) const;

        //- Inequality operator, compare %digest
        inline bool operator!=(const char* hexdigits) const;
};


// IOstream Operators

//- Output the %digest
inline Ostream& operator<<(Ostream& os, const SHA1& sha);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SHA1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
