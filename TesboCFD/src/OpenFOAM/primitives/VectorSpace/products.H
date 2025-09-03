/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Traits classes for inner and outer products of primitives.

\*---------------------------------------------------------------------------*/

#ifndef Foam_products_H
#define Foam_products_H

#include "direction.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Abstract template class to provide the form resulting from
//- the inner-product of two forms
template<class Cmpt, class Form1, class Form2>
class typeOfInnerProduct
{};

//- Abstract template class to provide the form resulting from
//- the outer-product of two forms
template<class Cmpt, class Form1, class Form2>
class typeOfOuterProduct
{};

//- Abstract template class to provide the transpose form of a form
template<class Cmpt, class Form>
class typeOfTranspose
{};


template<class Cmpt, direction rank>
class typeOfRank
{};


template<class Cmpt, direction rank>
class symmTypeOfRank
{};


//- The extended precision type (solveScalar for float)
template<class Type>
class typeOfSolve
{
public:

    typedef Type type;
};


//- The magnitude type for given argument.
template<class arg1>
class typeOfMag
{
public:

    typedef typename pTraits<typename pTraits<arg1>::cmptType>::magType type;
};


template<class arg1, class arg2>
class typeOfSum
{
public:

    typedef arg1 type;
};


template<class arg1, class arg2>
class outerProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg1>::cmptType,
        direction(pTraits<arg1>::rank) + direction(pTraits<arg2>::rank)
    >::type type;
};


//- Outer-product of identical types
template<class arg1>
class outerProduct1
:
    public outerProduct<arg1, arg1>
{
public:
};


template<class arg1, class arg2>
class crossProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg2>::cmptType,
        direction(pTraits<arg1>::rank) + direction(pTraits<arg2>::rank) - 1
    >::type type;
};

template<class arg1, class arg2>
class innerProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg1>::cmptType,
        direction(pTraits<arg1>::rank) + direction(pTraits<arg2>::rank) - 2
    >::type type;
};

template<class arg1, class arg2>
class scalarProduct
{
public:

    typedef typename pTraits<arg1>::cmptType type;
};


template<class arg1, direction arg2>
class powProduct
{
public:

    typedef typename symmTypeOfRank
    <
        typename pTraits<arg1>::cmptType,
        arg2*direction(pTraits<arg1>::rank)
    >::type type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Partial Specializations for non-VectorSpace quantities

template<class Cmpt>
class typeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
