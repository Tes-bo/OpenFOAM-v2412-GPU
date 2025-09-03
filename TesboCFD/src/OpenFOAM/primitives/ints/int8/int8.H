/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    int8_t

Description
    8bit signed integer

SourceFiles
    int8.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_int8_H
#define Foam_primitives_int8_H

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
inline word name(const int8_t val)
{
    return word(std::to_string(int(val)), false); // Needs no stripping
}


//- A word representation of uint8 value
template<>
struct nameOp<int8_t>
{
    word operator()(const int8_t val) const
    {
        return word(std::to_string(int(val)), false); // Needs no stripping
    }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read from stream (as label)
Istream& operator>>(Istream& is, int8_t& val);

//- Write to stream (as label)
Ostream& operator<<(Ostream& os, const int8_t val);


/*---------------------------------------------------------------------------*\
                        Specialization pTraits<int8_t>
\*---------------------------------------------------------------------------*/

//- Template specialization for pTraits<int8_t>
template<>
class pTraits<int8_t>
{
    int8_t p_;

public:

    // Typedefs

        //- Component type
        typedef int8_t cmptType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of int8_t is 0
        static constexpr direction rank = 0;

        //- Number of components in int8_t is 1
        static constexpr direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const int8_t zero;
        static const int8_t one;
        static const int8_t min;
        static const int8_t max;
        static const int8_t rootMax;
        static const int8_t rootMin;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(int8_t val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator int8_t() const noexcept { return p_; }

        //- Access the value
        operator int8_t&() noexcept { return p_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
