/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::expressions

Description
    Operations involving expressions.

SourceFiles
    exprOps.H

\*---------------------------------------------------------------------------*/

#ifndef exprOps_H
#define exprOps_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Namespace expressions Declarations
\*---------------------------------------------------------------------------*/

namespace expressions
{

//- Convert [0-1] values (usually scalars) as false/true
template<class T>
struct boolOp
{
    bool operator()(const T& val) const
    {
        return (0.5 < Foam::mag(val));
    }
};

//- Inverse of boolOp()
template<class T>
struct boolNotOp
{
    bool operator()(const T& val) const
    {
        return (Foam::mag(val) < 0.5);
    }
};

//- Convert to bool, pass-through for bool
template<>
struct boolOp<bool>
{
    bool operator()(bool val) const
    {
        return val;
    }
};

//- Inverse of boolOp, no conversion required for bool
template<>
struct boolNotOp<bool>
{
    bool operator()(bool val) const
    {
        return !val;
    }
};

//- Logical 'and', possibly with scalars as (false, true) values
template<class T>
struct boolAndOp
{
    bool operator()(const T& a, const T& b) const
    {
        return (boolOp<T>()(a) && boolOp<T>()(b));
    }
};

//- Logical 'or', possibly with scalars as (false, true) values
template<class T>
struct boolOrOp
{
    bool operator()(const T& a, const T& b) const
    {
        return (boolOp<T>()(a) || boolOp<T>()(b));
    }
};

//- Logical 'xor', possibly with scalars as (false, true) values
template<class T>
struct boolXorOp
{
    bool operator()(const T& a, const T& b) const
    {
        return (boolOp<T>()(a) != boolOp<T>()(b));
    }
};


} // End namespace expressions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
