/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Tout, class T1, class UnaryOp>
void Foam::FieldOps::assign
(
    Field<Tout>& result,
    const Field<T1>& a,
    const UnaryOp& op
)
{
    #ifdef FULLDEBUG
    if (result.size() != a.size())
    {
        FatalErrorInFunction
            << "Field sizes do not match: " << result.size() << " ("
            << a.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    std::transform(a.cbegin(), a.cend(), result.begin(), op);
}


template<class Tout, class T1, class T2, class BinaryOp>
void Foam::FieldOps::assign
(
    Field<Tout>& result,
    const Field<T1>& a,
    const Field<T2>& b,
    const BinaryOp& bop
)
{
    #ifdef FULLDEBUG
    if (result.size() != a.size() || result.size() != b.size())
    {
        FatalErrorInFunction
            << "Field sizes do not match: " << result.size() << " ("
            << a.size() << ' ' << b.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    std::transform(a.cbegin(), a.cend(), b.cbegin(), result.begin(), bop);
}


template<class T, class BinaryOp>
void Foam::FieldOps::ternary
(
    Field<T>& result,
    const Field<T>& a,
    const Field<T>& b,
    const BinaryOp& bop
)
{
    #ifdef FULLDEBUG
    if (result.size() != a.size() || result.size() != b.size())
    {
        FatalErrorInFunction
            << "Field sizes do not match: " << result.size() << " ("
            << a.size() << ' ' << b.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    forAll(result, i)
    {
        result[i] = bop(a[i], b[i]) ? a[i] : b[i];
    }
}


template<class T, class BoolListType, class FlipOp>
void Foam::FieldOps::ternarySelect
(
    Field<T>& result,
    const BoolListType& cond,
    const Field<T>& a,
    const Field<T>& b,
    const FlipOp& flip
)
{
    #ifdef FULLDEBUG
    if (result.size() != a.size() || result.size() != b.size())
    {
        FatalErrorInFunction
            << "Field sizes do not match: " << result.size() << " ("
            << a.size() << ' ' << b.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    forAll(result, i)
    {
        result[i] = flip(cond[i]) ? a[i] : b[i];
    }
}


template<class T, class FlipOp>
void Foam::FieldOps::ternarySelect
(
    Field<T>& result,
    const bitSet& cond,
    const Field<T>& a,
    const Field<T>& b,
    const FlipOp& flip
)
{
    #ifdef FULLDEBUG
    if (result.size() != a.size() || result.size() != b.size())
    {
        FatalErrorInFunction
            << "Field sizes do not match: " << result.size() << " ("
            << a.size() << ' ' << b.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    forAll(result, i)
    {
        result[i] = flip(cond[i]) ? a[i] : b[i];
    }
}


template<class T1, class T2>
Foam::Tuple2<T1,T2> Foam::FieldOps::findMinData
(
    const Field<T1>& vals,
    const Field<T2>& data
)
{
    Tuple2<T1,T2> result(pTraits<T1>::max, Zero);

    const label i = findMin(vals);
    if (i != -1)
    {
        result.first()  = vals[i];
        result.second() = data[i];
    }

    Pstream::combineReduce(result, minFirstEqOp<T1>());
    return result;
}


template<class T1, class T2>
Foam::Tuple2<T1,T2> Foam::FieldOps::findMaxData
(
    const Field<T1>& vals,
    const Field<T2>& data
)
{
    Tuple2<T1,T2> result(pTraits<T1>::min, Zero);

    const label i = findMax(vals);
    if (i != -1)
    {
        result.first()  = vals[i];
        result.second() = data[i];
    }

    Pstream::combineReduce(result, maxFirstEqOp<T1>());
    return result;
}


// ************************************************************************* //
