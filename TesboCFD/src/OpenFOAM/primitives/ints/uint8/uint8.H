/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    uint8_t

Description
    8bit unsigned integer

SourceFiles
    uint8.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_uint8_H
#define Foam_primitives_uint8_H

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

//- A word representation of uint8 value
inline word name(const uint8_t val)
{
    return word(std::to_string(int(val)), false); // Needs no stripping
}


//- A word representation of uint8 value
template<>
struct nameOp<uint8_t>
{
    word operator()(const uint8_t val) const
    {
        return word(std::to_string(int(val)), false); // Needs no stripping
    }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read uint8_t from stream.
uint8_t readUint8(Istream& is);

// IO operators are identical to direction, which is uint8_t


/*---------------------------------------------------------------------------*\
                        Specialization pTraits<uint8_t>
\*---------------------------------------------------------------------------*/

//- Template specialization for pTraits<uint8_t>
template<>
class pTraits<uint8_t>
{
    uint8_t p_;

public:

    // Typedefs

        //- Component type
        typedef uint8_t cmptType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of uint8_t is 0
        static constexpr direction rank = 0;

        //- Number of components in uint8_t is 1
        static constexpr direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const uint8_t zero;
        static const uint8_t one;
        static const uint8_t min;
        static const uint8_t max;
        static const uint8_t rootMax;
        static const uint8_t rootMin;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(uint8_t val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator uint8_t() const noexcept
        {
            return p_;
        }

        //- Access the value
        operator uint8_t&() noexcept
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
