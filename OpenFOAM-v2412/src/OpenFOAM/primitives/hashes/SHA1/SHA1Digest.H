/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::SHA1Digest

Description
    The SHA1 message digest.

See also
    Foam::SHA1

SourceFiles
    SHA1Digest.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SHA1Digest_H
#define Foam_SHA1Digest_H

#include <array>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class SHA1;

/*---------------------------------------------------------------------------*\
                         Class SHA1Digest Declaration
\*---------------------------------------------------------------------------*/

class SHA1Digest
{
    // Private Data

        //- The digest contents, which has 20 (uncoded) bytes
        std::array<unsigned char, 20> dig_;


    // Private Member Functions

        // Permit SHA1 to calculate the digest
        friend class SHA1;

        //- Raw digest data (20 bytes). Non-const access for SHA1
        unsigned char* data() noexcept { return dig_.data(); }


public:

    // Static Data Members

        //- A null digest (ie, all zero)
        static const SHA1Digest null;


    // Constructors

        //- Default construct a zero digest
        SHA1Digest();

        //- Read construct a digest from stringified content
        explicit SHA1Digest(Istream& is);

        //- Construct digest from raw or stringified content.
        //- The length is 20 for raw digest content and 40 (or 41) for
        //- stringified versions.
        SHA1Digest(const char* content, std::size_t len);

        //- Construct digest from raw or stringified content.
        //- The length is 20 for raw digest content and 40 (or 41) for
        //- stringified versions.
        SHA1Digest(const unsigned char* content, std::size_t len);


    // Member Functions

        //- Reset the digest to zero
        void clear();

        //- Return true if the digest is empty (ie, all zero).
        bool empty() const;

        //- Return (40-byte) text representation, optionally with '_' prefix
        std::string str(const bool prefixed=false) const;

        //- Read (40-byte) text representation.
        //  Since leading and intermediate underscores are skipped, a '_' can
        //  be prefixed to the text representation to use an unquoted
        //  SHA1Digest without parsing ambiguities as a number.
        Istream& read(Istream& is);

        //- Write (40-byte) text representation, optionally with '_' prefix
        Ostream& write(Ostream& os, const bool prefixed=false) const;


    // Low-level access

        //- Raw digest data (20 bytes) - const access
        const unsigned char* cdata() const noexcept { return dig_.data(); }

        //- Raw digest char data (20 bytes) - const access.
        //- For consistency with other objects, these are \em not unsigned.
        const char* cdata_bytes() const noexcept
        {
            return reinterpret_cast<const char*>(dig_.data());
        }

        //- Raw digest char data (20 bytes) - non-const access.
        //- For consistency with other objects, these are \em not unsigned.
        //- Use with caution - generally for broadcasting only.
        char* data_bytes() noexcept
        {
            return reinterpret_cast<char*>(dig_.data());
        }

        //- The number of bytes in digest (20)
        static constexpr unsigned size_bytes() noexcept { return 20; }

        //- The dimensioned size of the digest is always 20 bytes
        static constexpr unsigned max_size() noexcept { return 20; }


    // Member Operators

        //- Equality operator
        bool operator==(const SHA1Digest& rhs) const;

        //- Compare to (40-byte) text representation (eg, from sha1sum)
        //  An %empty string is equivalent to
        //  "0000000000000000000000000000000000000000"
        //  The hexdigits may optionally start with a '_' prefix
        bool operator==(const std::string& hexdigits) const;

        //- Compare to (40-byte) text representation (eg, from sha1sum)
        //  A %null or %empty string is equivalent to
        //  "0000000000000000000000000000000000000000"
        //  The hexdigits may optionally start with a '_' prefix
        bool operator==(const char* hexdigits) const;


        //- Inequality operator
        bool operator!=(const SHA1Digest& rhs) const;

        //- Inequality operator
        bool operator!=(const std::string& hexdigits) const;

        //- Inequality operator
        bool operator!=(const char* hexdigits) const;
};


// IOstream Operators

//- Read (40-byte) text representation (ignoring leading '_' prefix)
Istream& operator>>(Istream& is, SHA1Digest& dig);

//- Write (40-byte) text representation, unquoted and without prefix
Ostream& operator<<(Ostream& os, const SHA1Digest& dig);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
