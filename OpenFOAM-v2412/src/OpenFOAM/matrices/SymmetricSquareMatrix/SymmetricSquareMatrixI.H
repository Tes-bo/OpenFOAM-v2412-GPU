/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#define CHECK_MATRIX_IS_SQUARE(a, b)                                           \
    if (a != b)                                                                \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "Attempt to create a non-square matrix ("                       \
            << a << ", " << b << ')' << nl << abort(FatalError);               \
    }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix(const label n)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n)
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label n,
    const Foam::zero
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n, Foam::zero{})
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label n,
    const Type& val
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n, val)
{}


template<class Type>
template<class AnyType>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label n,
    const Identity<AnyType>
)
:
    SymmetricSquareMatrix<Type>(n, Foam::zero{})
{
    for (label i=0; i < n; ++i)
    {
        this->operator()(i, i) = pTraits<Type>::one;
    }
}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix(Istream& is)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(is)
{
    CHECK_MATRIX_IS_SQUARE(this->m(), this->n());
}


template<class Type>
inline Foam::autoPtr<Foam::SymmetricSquareMatrix<Type>>
Foam::SymmetricSquareMatrix<Type>::clone() const
{
    return autoPtr<SymmetricSquareMatrix<Type>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::SymmetricSquareMatrix<Type>::operator=(const Foam::zero)
{
    Matrix<SymmetricSquareMatrix<Type>, Type>::operator=(Foam::zero{});
}


template<class Type>
inline void Foam::SymmetricSquareMatrix<Type>::operator=(const Type& val)
{
    Matrix<SymmetricSquareMatrix<Type>, Type>::operator=(val);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef CHECK_MATRIX_IS_SQUARE

// ************************************************************************* //
