/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "DiagonalMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label n)
:
    List<Type>(n)
{}


template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label n, const Foam::zero)
:
    List<Type>(n, Foam::zero{})
{}


template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label n, const Type& val)
:
    List<Type>(n, val)
{}


template<class Type>
template<class Form>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const Matrix<Form, Type>& mat)
:
    List<Type>(min(mat.m(), mat.n()))
{
    label i = 0;

    for (Type& val : *this)
    {
        val = mat(i, i);
        ++i;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::DiagonalMatrix<Type>::invert()
{
    for (Type& val : *this)
    {
        // If mag(val)<VSMALL, leave untouched
        // forcing of val=Zero sometimes confuses compiler
        if (mag(val) > VSMALL)
        {
            val = Type(1)/val;
        }
    }
}


template<class Type>
template<class CompOp>
Foam::List<Foam::label> Foam::DiagonalMatrix<Type>::sortPermutation
(
    CompOp& compare
) const
{
    List<label> p(this->size());
    std::iota(p.begin(), p.end(), 0);
    std::sort
    (
        p.begin(),
        p.end(),
        [&](label i, label j){ return compare((*this)[i], (*this)[j]); }
    );

    return p;
}


template<class Type>
void Foam::DiagonalMatrix<Type>::applyPermutation(const List<label>& p)
{
    #ifdef FULLDEBUG
    if (this->size() != p.size())
    {
        FatalErrorInFunction
            << "Attempt to column-reorder according to an uneven list: " << nl
            << "DiagonalMatrix diagonal size = " << this->size() << nl
            << "Permutation list size = " << p.size() << nl
            << abort(FatalError);
    }
    #endif

    List<bool> pass(p.size(), false);

    for (label i = 0; i < p.size(); ++i)
    {
        if (pass[i])
        {
            continue;
        }
        pass[i] = true;
        label prev = i;
        label j = p[i];
        while (i != j)
        {
            Foam::Swap((*this)[prev], (*this)[j]);
            pass[j] = true;
            prev = j;
            j = p[j];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the matrix inverse as a DiagonalMatrix
template<class Type>
DiagonalMatrix<Type> inv(const DiagonalMatrix<Type>& mat)
{
    // Construct with fall-back value conditional calculation
    // of invert to avoid over-eager compiler optimisation
    DiagonalMatrix<Type> Ainv(mat.size(), Zero);

    Type* iter = Ainv.begin();

    for (const Type& val : mat)
    {
        if (mag(val) > VSMALL)
        {
            *iter = Type(1)/val;
        }

        ++iter;
    }

    return Ainv;
}


//- Return Matrix column-reordered according to
//- a given permutation labelList
template<class Type>
DiagonalMatrix<Type> applyPermutation
(
    const DiagonalMatrix<Type>& mat,
    const List<label>& p
)
{
    #ifdef FULLDEBUG
    if (mat.size() != p.size())
    {
        FatalErrorInFunction
            << "Attempt to column-reorder according to an uneven list: " << nl
            << "DiagonalMatrix diagonal size = " << mat.size() << nl
            << "Permutation list size = " << p.size() << nl
            << abort(FatalError);
    }
    #endif

    DiagonalMatrix<Type> reordered(mat.size());

    label j = 0;
    for (const label i : p)
    {
        reordered[j] = mat[i];
        ++j;
    }

    return reordered;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
