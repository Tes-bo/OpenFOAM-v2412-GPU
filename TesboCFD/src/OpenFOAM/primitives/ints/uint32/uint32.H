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
    uint32_t

Description
    32bit unsigned integer

SourceFiles
    uint32.C
    uint32IO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_uint32_H
#define Foam_primitives_uint32_H

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

//- A word representation of uint32 value
inline word name(const uint32_t val)
{
    return word(std::to_string(val), false); // Needs no stripping
}

//- A word representation of uint32 value
template<>
struct nameOp<uint32_t>
{
    word operator()(const uint32_t val) const
    {
        return word(std::to_string(val), false); // Needs no stripping
    }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read uint32_t from stream
uint32_t readUint32(Istream& is);

//- Parse entire buffer as a uint32_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
uint32_t readUint32(const char* buf);

//- Parse entire string as a uint32_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline uint32_t readUint32(const std::string& str)
{
    return readUint32(str.c_str());
}

//- Read entire buffer as a uint32_t, skipping leading/trailing whitespace.
//  \return True if successful.
bool readUint32(const char* buf, uint32_t& val);

//- Read entire string as a uint32_t, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool readUint32(const std::string& str, uint32_t& val)
{
    return readUint32(str.c_str(), val);
}

//- Same as readUint32
//  \return True if successful.
inline bool read(const char* buf, uint32_t& val)
{
    return readUint32(buf, val);
}

//- Same as readUint32
//  \return True if successful.
inline bool read(const std::string& str, uint32_t& val)
{
    return readUint32(str, val);
}


Istream& operator>>(Istream& is, uint32_t& val);
Ostream& operator<<(Ostream& os, const uint32_t val);


/*---------------------------------------------------------------------------*\
                       Specialization pTraits<uint32_t>
\*---------------------------------------------------------------------------*/

//- Template specialization for pTraits<uint32_t>
template<>
class pTraits<uint32_t>
{
    uint32_t p_;

public:

    // Typedefs

        //- Component type
        typedef uint32_t cmptType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of uint32_t is 0
        static constexpr direction rank = 0;

        //- Number of components in uint32_t is 1
        static constexpr direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const uint32_t zero;
        static const uint32_t one;
        static const uint32_t min;
        static const uint32_t max;
        static const uint32_t rootMax;
        static const uint32_t rootMin;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(uint32_t val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator uint32_t() const noexcept { return p_; }

        //- Access the value
        operator uint32_t&() noexcept { return p_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
