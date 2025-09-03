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

#define CHECK_MATRIX_IS_SQUARE(a, b)                                           \
    if (a != b)                                                                \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "Attempt to create a non-square matrix ("                       \
            << a << ", " << b << ')' << nl << abort(FatalError);               \
    }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix(const label n)
:
    Matrix<SquareMatrix<Type>, Type>(n, n)
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label n,
    const Foam::zero
)
:
    Matrix<SquareMatrix<Type>, Type>(n, n, Foam::zero{})
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label n,
    const Type& val
)
:
    Matrix<SquareMatrix<Type>, Type>(n, n, val)
{}


template<class Type>
template<class AnyType>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label n,
    const Identity<AnyType>
)
:
    SquareMatrix<Type>(n, Foam::zero{})
{
    for (label i = 0; i < n; ++i)
    {
        this->operator()(i, i) = pTraits<Type>::one;
    }
}


template<class Type>
template<class AnyType>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const labelPair& dims,
    const Identity<AnyType>
)
:
    Matrix<SquareMatrix<Type>, Type>(dims, Foam::zero{})
{
    CHECK_MATRIX_IS_SQUARE(dims.first(), dims.second());

    for (label i = 0; i < dims.first(); ++i)
    {
        this->operator()(i, i) = pTraits<Type>::one;
    }
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const labelPair& dims
)
:
    Matrix<SquareMatrix<Type>, Type>(dims)
{
    CHECK_MATRIX_IS_SQUARE(dims.first(), dims.second());
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const labelPair& dims,
    const Foam::zero
)
:
    Matrix<SquareMatrix<Type>, Type>(dims, Foam::zero{})
{
    CHECK_MATRIX_IS_SQUARE(dims.first(), dims.second());
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const labelPair& dims,
    const Type& val
)
:
    Matrix<SquareMatrix<Type>, Type>(dims, val)
{
    CHECK_MATRIX_IS_SQUARE(dims.first(), dims.second());
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label m,
    const label n,
    const Foam::zero
)
:
    Matrix<SquareMatrix<Type>, Type>(m, n, Foam::zero{})
{
    CHECK_MATRIX_IS_SQUARE(m, n);
}


template<class Type>
template<class MatrixType>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const ConstMatrixBlock<MatrixType>& mat
)
:
    Matrix<SquareMatrix<Type>, Type>(mat)
{
    // Check is square?
}


template<class Type>
template<class MatrixType>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const MatrixBlock<MatrixType>& mat
)
:
    Matrix<SquareMatrix<Type>, Type>(mat)
{
    // Check is square?
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const RectangularMatrix<Type>& mat
)
:
    Matrix<SquareMatrix<Type>, Type>(mat)
{
    CHECK_MATRIX_IS_SQUARE(mat.m(), mat.n());
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix(Istream& is)
:
    Matrix<SquareMatrix<Type>, Type>(is)
{
    CHECK_MATRIX_IS_SQUARE(this->m(), this->n());
}


template<class Type>
inline Foam::autoPtr<Foam::SquareMatrix<Type>>
Foam::SquareMatrix<Type>::clone() const
{
    return autoPtr<SquareMatrix<Type>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::SquareMatrix<Type>::resize(const label m)
{
    Matrix<SquareMatrix<Type>, Type>::resize(m, m);
}


template<class Type>
inline void Foam::SquareMatrix<Type>::resize_nocopy(const label m)
{
    Matrix<SquareMatrix<Type>, Type>::resize_nocopy(m, m);
}


template<class Type>
inline void Foam::SquareMatrix<Type>::resize(const label m, const label n)
{
    if (m != n)
    {
        FatalErrorInFunction<< "Disallowed use of resize() for SquareMatrix"
            << abort(FatalError);
    }

    Matrix<SquareMatrix<Type>, Type>::resize(m, m);
}


template<class Type>
inline void Foam::SquareMatrix<Type>::setSize(const label m)
{
    Matrix<SquareMatrix<Type>, Type>::resize(m, m);
}


template<class Type>
inline void Foam::SquareMatrix<Type>::shallowResize(const label m)
{
    Matrix<SquareMatrix<Type>, Type>::shallowResize(m, m);
}


template<class Type>
inline bool Foam::SquareMatrix<Type>::symmetric() const
{
    for (label n = 0; n < this->n() - 1; ++n)
    {
        for (label m = this->m() - 1; n < m; --m)
        {
            if (SMALL < mag((*this)(n, m) - (*this)(m, n)))
            {
                return false;
            }
        }
    }
    return true;
}


template<class Type>
inline bool Foam::SquareMatrix<Type>::tridiagonal() const
{
    for (label i = 0; i < this->m(); ++i)
    {
        for (label j = 0; j < this->n(); ++j)
        {
            const Type& val = (*this)(i, j);

            if ((i == j) || (i - 1 == j) || (i + 1 == j))
            {
                if (mag(val) < SMALL)
                {
                    return false;
                }
            }
            else if (SMALL < mag(val))
            {
                return false;
            }
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::SquareMatrix<Type>::operator=(SquareMatrix<Type>&& mat)
{
    this->transfer(mat);
}


template<class Type>
inline void Foam::SquareMatrix<Type>::operator=(const Foam::zero)
{
    Matrix<SquareMatrix<Type>, Type>::operator=(Foam::zero{});
}


template<class Type>
inline void Foam::SquareMatrix<Type>::operator=(const Type& val)
{
    Matrix<SquareMatrix<Type>, Type>::operator=(val);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// Return the outer product of Field-Field as SquareMatrix
template<class Type>
inline Foam::SquareMatrix<Type> symmOuter
(
    const Field<Type>& f1,
    const Field<Type>& f2
)
{
    SquareMatrix<Type> f1f2T(f1.size());

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef CHECK_MATRIX_IS_SQUARE

// ************************************************************************* //
