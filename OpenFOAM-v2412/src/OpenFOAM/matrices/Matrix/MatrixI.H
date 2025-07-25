/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "MatrixBlock.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::doAlloc()
{
    const label len = size();

    if (len > 0)
    {
        // With sign-check to avoid spurious -Walloc-size-larger-than
        v_ = new Type[len];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix() noexcept
:
    mRows_(0),
    nCols_(0),
    v_(nullptr)
{}


template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix(const labelPair& dims)
:
    Matrix<Form, Type>(dims.first(), dims.second())
{}


template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix(const labelPair& dims, const Foam::zero)
:
    Matrix<Form, Type>(dims.first(), dims.second(), Foam::zero{})
{}


template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix(const labelPair& dims, const Type& val)
:
    Matrix<Form, Type>(dims.first(), dims.second(), val)
{}


template<class Form, class Type>
inline Foam::autoPtr<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::clone() const
{
    return autoPtr<Matrix<Form, Type>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::size() const
{
    return mRows_ * nCols_;
}


template<class Form, class Type>
inline Foam::labelPair Foam::Matrix<Form, Type>::sizes() const
{
    return labelPair(mRows_, nCols_);
}


template<class Form, class Type>
inline bool Foam::Matrix<Form, Type>::empty() const noexcept
{
    return !mRows_ || !nCols_;
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checki(const label i) const
{
    if (!mRows_ || !nCols_)
    {
        FatalErrorInFunction
            << "Attempt to access element from empty matrix"
            << abort(FatalError);
    }
    if (i < 0 || mRows_ <= i)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << mRows_-1
            << abort(FatalError);
    }
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checkj(const label j) const
{
    if (!mRows_ || !nCols_)
    {
        FatalErrorInFunction
            << "Attempt to access element from empty matrix"
            << abort(FatalError);
    }
    if (j < 0 || nCols_ <= j)
    {
        FatalErrorInFunction
            << "index " << j << " out of range 0 ... " << nCols_-1
            << abort(FatalError);
    }
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checkSize() const
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect size (" << mRows_ << ", " << nCols_ << ')' << nl
            << abort(FatalError);
    }
    // Could also check for odd sizes, like (0, N) and make (0,0)
}


template<class Form, class Type>
inline bool Foam::Matrix<Form, Type>::uniform() const
{
    const label len = size();

    if (!len)
    {
        return false;
    }

    // std::all_of()
    for (label idx = 1; idx < len; ++idx)
    {
        if (v_[0] != v_[idx])
        {
            return false;
        }
    }

    return true;
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::cdata() const noexcept
{
    return v_;
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::data() noexcept
{
    return v_;
}


template<class Form, class Type>
inline const char* Foam::Matrix<Form, Type>::cdata_bytes() const noexcept
{
    return reinterpret_cast<const char*>(v_);
}


template<class Form, class Type>
inline char* Foam::Matrix<Form, Type>::data_bytes() noexcept
{
    return reinterpret_cast<char*>(v_);
}


template<class Form, class Type>
inline std::streamsize Foam::Matrix<Form, Type>::size_bytes() const noexcept
{
    return mRows_*nCols_*sizeof(Type);
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::rowData(const label irow) const
{
    #ifdef FULLDEBUG
    checki(irow);
    #endif
    return (v_ + irow*nCols_);
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::rowData(const label irow)
{
    #ifdef FULLDEBUG
    checki(irow);
    #endif
    return (v_ + irow*nCols_);
}


template<class Form, class Type>
inline const Type& Foam::Matrix<Form, Type>::at(const label idx) const
{
    #ifdef FULLDEBUG
    if (idx < 0 || this->size() <= idx)
    {
        FatalErrorInFunction
            << "index " << idx << " out of range 0 ... " << this->size()
            << abort(FatalError);
    }
    #endif
    return *(v_ + idx);
}


template<class Form, class Type>
inline Type& Foam::Matrix<Form, Type>::at(const label idx)
{
    #ifdef FULLDEBUG
    if (idx < 0 || this->size() <= idx)
    {
        FatalErrorInFunction
            << "index " << idx << " out of range 0 ... " << this->size()
            << abort(FatalError);
    }
    #endif
    return *(v_ + idx);
}


template<class Form, class Type>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::subColumn
(
    const label colIndex,
    const label rowIndex,
    label len
) const
{
    if (len < 0)
    {
        len = mRows_ - rowIndex;
    }

    return ConstMatrixBlock<mType>
    (
        *this,
        len, // rows
        1,
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::subRow
(
    const label rowIndex,
    const label colIndex,
    label len
) const
{
    if (len < 0)
    {
        len = nCols_ - colIndex;
    }

    return ConstMatrixBlock<mType>
    (
        *this,
        1,
        len, // columns
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::subMatrix
(
    const label rowIndex,
    const label colIndex,
    label szRows,
    label szCols
) const
{
    if (szRows < 0) szRows = mRows_ - rowIndex;
    if (szCols < 0) szCols = nCols_ - colIndex;

    return ConstMatrixBlock<mType>
    (
        *this,
        szRows,
        szCols,
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
template<class VectorSpace>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::block
(
    const label rowIndex,
    const label colIndex
) const
{
    return ConstMatrixBlock<mType>
    (
        *this,
        VectorSpace::mRows,
        VectorSpace::nCols,
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::subColumn
(
    const label colIndex,
    const label rowIndex,
    label len
)
{
    if (len < 0)
    {
        len = mRows_ - rowIndex;
    }

    return MatrixBlock<mType>
    (
        *this,
        len, // rows
        1,
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::subRow
(
    const label rowIndex,
    const label colIndex,
    label len
)
{
    if (len < 0)
    {
        len = nCols_ - colIndex;
    }

    return MatrixBlock<mType>
    (
        *this,
        1,
        len, // columns
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::subMatrix
(
    const label rowIndex,
    const label colIndex,
    label szRows,
    label szCols
)
{
    if (szRows < 0) szRows = mRows_ - rowIndex;
    if (szCols < 0) szCols = nCols_ - colIndex;

    return MatrixBlock<mType>
    (
        *this,
        szRows,
        szCols,
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
template<class VectorSpace>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::block
(
    const label rowIndex,
    const label colIndex
)
{
    return MatrixBlock<mType>
    (
        *this,
        VectorSpace::mRows,
        VectorSpace::nCols,
        rowIndex,
        colIndex
    );
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::setSize(const label m, const label n)
{
    resize(m, n);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::shallowResize(const label m, const label n)
{
    mRows_ = m;
    nCols_ = n;
}


template<class Form, class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::Amul
(
    const UList<Type>& x
) const
{
    return this->AmulImpl(x);
}


template<class Form, class Type>
template<class Addr>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::Amul
(
    const IndirectListBase<Type, Addr>& x
) const
{
    return this->AmulImpl(x);
}


template<class Form, class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::Tmul
(
    const UList<Type>& x
) const
{
    return this->TmulImpl(x);
}


template<class Form, class Type>
template<class Addr>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::Tmul
(
    const IndirectListBase<Type, Addr>& x
) const
{
    return this->TmulImpl(x);
}


// * * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * //

template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::iterator
Foam::Matrix<Form, Type>::begin() noexcept
{
    return v_;
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::iterator
Foam::Matrix<Form, Type>::end() noexcept
{
    return v_ + (mRows_ * nCols_);
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::cbegin() const noexcept
{
    return v_;
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::cend() const noexcept
{
    return v_ + (mRows_ * nCols_);
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::begin() const noexcept
{
    return v_;
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::end() const noexcept
{
    return v_ + (mRows_ * nCols_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
inline const Type& Foam::Matrix<Form, Type>::operator()
(
    const label irow,
    const label jcol
) const
{
    #ifdef FULLDEBUG
    checki(irow);
    checkj(jcol);
    #endif
    return v_[irow*nCols_ + jcol];
}


template<class Form, class Type>
inline Type& Foam::Matrix<Form, Type>::operator()
(
    const label irow,
    const label jcol
)
{
    #ifdef FULLDEBUG
    checki(irow);
    checkj(jcol);
    #endif
    return v_[irow*nCols_ + jcol];
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::operator[](const label irow) const
{
    #ifdef FULLDEBUG
    checki(irow);
    #endif
    return v_ + irow*nCols_;
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::operator[](const label irow)
{
    #ifdef FULLDEBUG
    checki(irow);
    #endif
    return v_ + irow*nCols_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Matrix-vector multiplication (A * x), where x is a column vector
template<class Form, class Type>
inline tmp<Field<Type>> operator*
(
    const Matrix<Form, Type>& mat,
    const UList<Type>& x
)
{
    return mat.Amul(x);
}


//- Matrix-vector multiplication (A * x), where x is a column vector
template<class Form, class Type, class Addr>
inline tmp<Field<Type>> operator*
(
    const Matrix<Form, Type>& mat,
    const IndirectListBase<Type, Addr>& x
)
{
    return mat.Amul(x);
}


//- Vector-Matrix multiplication (x * A), where x is a row vector
template<class Form, class Type>
inline tmp<Field<Type>> operator*
(
    const UList<Type>& x,
    const Matrix<Form, Type>& mat
)
{
    return mat.Tmul(x);
}


//- Vector-Matrix multiplication (x * A), where x is a row vector
template<class Form, class Type, class Addr>
inline tmp<Field<Type>> operator*
(
    const IndirectListBase<Type, Addr>& x,
    const Matrix<Form, Type>& mat
)
{
    return mat.Tmul(x);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
