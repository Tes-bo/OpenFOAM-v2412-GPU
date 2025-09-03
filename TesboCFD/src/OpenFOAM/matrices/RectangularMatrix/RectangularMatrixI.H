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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const label m,
    const label n
)
:
    Matrix<RectangularMatrix<Type>, Type>(m, n)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const label n
)
:
    Matrix<RectangularMatrix<Type>, Type>(n, n)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const label m,
    const label n,
    const Foam::zero
)
:
    Matrix<RectangularMatrix<Type>, Type>(m, n, Foam::zero{})
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const label m,
    const label n,
    const Type& val
)
:
    Matrix<RectangularMatrix<Type>, Type>(m, n, val)
{}


template<class Type>
template<class AnyType>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const labelPair& dims,
    const Identity<AnyType>
)
:
    RectangularMatrix<Type>(dims.first(), dims.second(), Foam::zero{})
{
    for (label i = 0; i < min(dims.first(), dims.second()); ++i)
    {
        this->operator()(i, i) = pTraits<Type>::one;
    }
}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const labelPair& dims
)
:
    RectangularMatrix<Type>(dims.first(), dims.second())
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const labelPair& dims,
    const Foam::zero
)
:
    RectangularMatrix<Type>(dims.first(), dims.second(), Foam::zero{})
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const labelPair& dims,
    const Type& val
)
:
    RectangularMatrix<Type>(dims.first(), dims.second(), val)
{}


template<class Type>
template<class MatrixType>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const ConstMatrixBlock<MatrixType>& mat
)
:
    Matrix<RectangularMatrix<Type>, Type>(mat)
{}


template<class Type>
template<class MatrixType>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const MatrixBlock<MatrixType>& mat
)
:
    Matrix<RectangularMatrix<Type>, Type>(mat)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const SquareMatrix<Type>& mat
)
:
    Matrix<RectangularMatrix<Type>, Type>(mat)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix(Istream& is)
:
    Matrix<RectangularMatrix<Type>, Type>(is)
{}


template<class Type>
inline Foam::autoPtr<Foam::RectangularMatrix<Type>>
Foam::RectangularMatrix<Type>::clone() const
{
    return autoPtr<RectangularMatrix<Type>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::RectangularMatrix<Type>::operator=
(
    RectangularMatrix<Type>&& mat
)
{
    this->transfer(mat);
}


template<class Type>
inline void Foam::RectangularMatrix<Type>::operator=(const Foam::zero)
{
    Matrix<RectangularMatrix<Type>, Type>::operator=(Foam::zero{});
}


template<class Type>
inline void Foam::RectangularMatrix<Type>::operator=(const Type& val)
{
    Matrix<RectangularMatrix<Type>, Type>::operator=(val);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
class typeOfInnerProduct<Type, RectangularMatrix<Type>, RectangularMatrix<Type>>
{
public:

    typedef RectangularMatrix<Type> type;
};


template<class Type>
class typeOfInnerProduct<Type, RectangularMatrix<Type>, SquareMatrix<Type>>
{
public:

    typedef RectangularMatrix<Type> type;
};


template<class Type>
class typeOfInnerProduct<Type, SquareMatrix<Type>, RectangularMatrix<Type>>
{
public:

    typedef RectangularMatrix<Type> type;
};


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// Return the outer product of Field-Field as RectangularMatrix
template<class Type>
inline Foam::RectangularMatrix<Type> outer
(
    const Field<Type>& f1,
    const Field<Type>& f2
)
{
    RectangularMatrix<Type> f1f2T(f1.size(), f2.size());

    for (label i = 0; i < f1f2T.m(); ++i)
    {
        for (label j = 0; j < f1f2T.n(); ++j)
        {
            f1f2T(i, j) = f1[i]*f2[j];
        }
    }

    return f1f2T;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
