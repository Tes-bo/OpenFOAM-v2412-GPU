/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Namespace
    Foam::Detail::ListPolicy

Description
    Additional compile-time controls of List behaviour

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListPolicy_H
#define Foam_ListPolicy_H

#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class keyType;
class word;
class wordRe;

namespace Detail
{
namespace ListPolicy
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Number of items before requiring line-breaks in the list output.
//
//  Default definition: 10
template<class T>
struct short_length : std::integral_constant<int,10> {};

// Can override on a per-type basis
// Eg,
// template<> struct short_length<label> : std::integral_constant<int,20> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Can suppress additional line breaks separate ASCII data content
//- when the data elements are primitives, or contiguous
//
//  Default definition: (integral | floating-point) are contiguous and thus
//  never need any line breaks
template<class T>
struct no_linebreak : std::is_arithmetic<T> {};

// Specialization for word-like classes
// These elements are normally fairly short, so ok to output a few (eg, 10)
// of them on a single line.

//- Suppress line-breaks for keyType
template<> struct no_linebreak<keyType> : std::true_type {};

//- Suppress line-breaks for word
template<> struct no_linebreak<word> : std::true_type {};

//- Suppress line-breaks for wordRe
template<> struct no_linebreak<wordRe> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Classification of list/container uniformity.
//- The values can be used with bit-wise \c or reduction
enum uniformity : unsigned char
{
    EMPTY = 0,          //!< An empty container
    UNIFORM = 0x1,      //!< Container (non-empty) with identical values
    NONUNIFORM = 0x2,   //!< Container (non-empty) with different values
    MIXED = 0x3         //!< Mixed uniform/non-uniform (eg, after reduction)
};

//- Algorithm to determine list/container uniformity
template<class InputIt>
enum uniformity check_uniformity(InputIt first, InputIt last)
{
    if (first == last) return uniformity::EMPTY;

    // Like std::all_of() with checking against element 0,
    // but without using a lambda with auto type (pre C++14) etc.

    const auto& elem0 = *first;

    for ((void)++first; (first != last); (void)++first)
    {
        if (elem0 != *first)
        {
            return uniformity::NONUNIFORM;
        }
    }

    return uniformity::UNIFORM;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ListPolicy
} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
