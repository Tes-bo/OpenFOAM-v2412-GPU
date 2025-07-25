/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "LLTMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LLTMatrix<Type>::LLTMatrix()
{}


template<class Type>
Foam::LLTMatrix<Type>::LLTMatrix(const SquareMatrix<Type>& mat)
{
    decompose(mat);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::LLTMatrix<Type>::decompose(const SquareMatrix<Type>& mat)
{
    SquareMatrix<Type>& LLT = *this;

    // Initialize the LLT decomposition matrix to M
    LLT = mat;

    const label m = LLT.m();

    for (label i = 0; i < m; ++i)
    {
        for (label j = 0; j < m; ++j)
        {
            if (j > i)
            {
                LLT(i, j) = Zero;
                continue;
            }

            Type sum = LLT(i, j);

            for (label k = 0; k < j; ++k)
            {
                sum -= LLT(i, k)*LLT(j, k);
            }

            if (i > j)
            {
                LLT(i, j) = sum/LLT(j, j);
            }
            else if (sum > 0)
            {
                LLT(i, i) = sqrt(sum);
            }
            else
            {
                FatalErrorInFunction
                    << "Cholesky decomposition failed, "
                       "matrix is not symmetric positive definite"
                    << abort(FatalError);
            }
        }
    }
}


template<class Type>
template<template<typename> class ListContainer>
void Foam::LLTMatrix<Type>::solveImpl
(
    List<Type>& x,
    const ListContainer<Type>& source
) const
{
    // If x and source are different, copy initialize x = source
    if (&x != &source)
    {
        x = source;
    }

    const SquareMatrix<Type>& LLT = *this;
    const label m = LLT.m();

    for (label i = 0; i < m; ++i)
    {
        Type sum = source[i];

        for (label j = 0; j < i; ++j)
        {
            sum = sum - LLT(i, j)*x[j];
        }

        x[i] = sum/LLT(i, i);
    }

    for (label i = m - 1; i >= 0; --i)
    {
        Type sum = x[i];

        for (label j = i + 1; j < m; ++j)
        {
            sum = sum - LLT(j, i)*x[j];
        }

        x[i] = sum/LLT(i, i);
    }
}


template<class Type>
void Foam::LLTMatrix<Type>::solve
(
    List<Type>& x,
    const UList<Type>& source
) const
{
    solveImpl(x, source);
}


template<class Type>
template<class Addr>
void Foam::LLTMatrix<Type>::solve
(
    List<Type>& x,
    const IndirectListBase<Type, Addr>& source
) const
{
    solveImpl(x, source);
}

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::LLTMatrix<Type>::solve
(
    const UList<Type>& source
) const
{
    auto tresult(tmp<Field<Type>>::New(source.size()));

    solve(tresult.ref(), source);

    return tresult;
}


template<class Type>
template<class Addr>
Foam::tmp<Foam::Field<Type>> Foam::LLTMatrix<Type>::solve
(
    const IndirectListBase<Type, Addr>& source
) const
{
    auto tresult(tmp<Field<Type>>::New(source.size()));

    solve(tresult.ref(), source);

    return tresult;
}


// ************************************************************************* //
