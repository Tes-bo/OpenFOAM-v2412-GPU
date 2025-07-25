/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

Primitive
    uint64_t

Description
    64bit unsigned integer

SourceFiles
    uint64.C
    uint64IO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_uint64_H
#define Foam_primitives_uint64_H

#include <cstdint>
#include <climits>
#include <cstdlib>

#include "direction.H"
#include "pTraits.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A word representation of uint64 value
inline word name(const uint64_t val)
{
    return word(std::to_string(val), false); // Needs no stripping
}


//- A word representation of uint64 value
template<>
struct nameOp<uint64_t>
{
    word operator()(const uint64_t val) const
    {
        return word(std::to_string(val), false); // Needs no stripping
    }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read uint64_t from stream.
uint64_t readUint64(Istream& is);

//- Parse entire buffer as uint64_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
uint64_t readUint64(const char* buf);

//- Parse entire string as uint64_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline uint64_t readUint64(const std::string& str)
{
    return readUint64(str.c_str());
}

//- Parse entire buffer as uint64_t, skipping leading/trailing whitespace.
//  \return True if successful.
bool readUint64(const char* buf, uint64_t& val);

//- Parse entire string as uint64_t, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool readUint64(const std::string& str, uint64_t& val)
{
    return readUint64(str.c_str(), val);
}

//- Same as readUint64
//  \return True if successful.
inline bool read(const char* buf, uint64_t& val)
{
    return readUint64(buf, val);
}

//- Same as readUint64
//  \return True if successful.
inline bool read(const std::string& str, uint64_t& val)
{
    return readUint64(str, val);
}


Istream& operator>>(Istream& is, uint64_t& val);
Ostream& operator<<(Ostream& os, const uint64_t val);

// On Darwin:
// unsigned long is not unambiguously (uint32_t | uint64_t)
// - explicitly resolve for output
#ifdef __APPLE__
    Ostream& operator<<(Ostream& os, const unsigned long val);
#endif


/*---------------------------------------------------------------------------*\
                       Specialization pTraits<uint64_t>
\*---------------------------------------------------------------------------*/

//- Template specialization for pTraits<uint64_t>
template<>
class pTraits<uint64_t>
{
    uint64_t p_;

public:

    // Typedefs

        //- Component type
        typedef uint64_t cmptType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of uint64_t is 0
        static constexpr direction rank = 0;

        //- Number of components in uint64_t is 1
        static constexpr direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const uint64_t zero;
        static const uint64_t one;
        static const uint64_t min;
        static const uint64_t max;
        static const uint64_t rootMax;
        static const uint64_t rootMin;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(uint64_t val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator uint64_t() const noexcept { return p_; }

        //- Access the value
        operator uint64_t&() noexcept { return p_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
