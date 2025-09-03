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
    int32_t

Description
    32bit signed integer

SourceFiles
    int32.C
    int32IO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_int32_H
#define Foam_primitives_int32_H

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

//- A word representation of int32 value
inline word name(const int32_t val)
{
    return word(std::to_string(val), false); // Needs no stripping
}


//- A word representation of int32 value
template<>
struct nameOp<int32_t>
{
    word operator()(const int32_t val) const
    {
        return word(std::to_string(val), false); // Needs no stripping
    }
};


inline int32_t mag(const int32_t val)
{
    return ::abs(val);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read int32_t from stream
int32_t readInt32(Istream& is);

//- Parse entire buffer as a int32_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
int32_t readInt32(const char* buf);

//- Parse entire string as a int32_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline int32_t readInt32(const std::string& str)
{
    return readInt32(str.c_str());
}

//- Read entire buffer as a int32_t, skipping leading/trailing whitespace.
//  \return True if successful.
bool readInt32(const char* buf, int32_t& val);

//- Read entire string as a int32_t, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool readInt32(const std::string& str, int32_t& val)
{
    return readInt32(str.c_str(), val);
}

//- Same as readInt32
//  \return True if successful.
inline bool read(const char* buf, int32_t& val)
{
    return readInt32(buf, val);
}

//- Same as readInt32
//  \return True if successful.
inline bool read(const std::string& str, int32_t& val)
{
    return readInt32(str, val);
}


Istream& operator>>(Istream& is, int32_t& val);
Ostream& operator<<(Ostream& os, const int32_t val);

// 32bit compilation with long as int32_t
// - resolve explicitly for input and output
//
// Test works for gcc, icc, llvm.
#if (__SIZEOF_LONG__ == 4)
    Istream& operator>>(Istream& is, long& val);
    Ostream& operator<<(Ostream& os, const long val);
#endif


/*---------------------------------------------------------------------------*\
                        Specialization pTraits<int32_t>
\*---------------------------------------------------------------------------*/

//- Template specialization for pTraits<int32_t>
template<>
class pTraits<int32_t>
{
    int32_t p_;

public:

    // Typedefs

        //- Component type
        typedef int32_t cmptType;

        //- Magnitude type
        typedef int32_t magType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of int32_t is 0
        static constexpr direction rank = 0;

        //- Number of components in int32_t is 1
        static constexpr direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const int32_t zero;
        static const int32_t one;
        static const int32_t min;
        static const int32_t max;
        static const int32_t rootMax;
        static const int32_t rootMin;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(int32_t val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator int32_t() const noexcept { return p_; }

        //- Access the value
        operator int32_t&() noexcept { return p_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
