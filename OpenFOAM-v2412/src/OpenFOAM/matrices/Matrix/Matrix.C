/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "Matrix.H"
#include <functional>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Form, class Type>
template<class ListType>
Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::AmulImpl
(
    const ListType& x
) const
{
    const Matrix<Form, Type>& mat = *this;

    #ifdef FULLDEBUG
    if (mat.n() != x.size())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible Matrix and Vector:" << nl
            << "Matrix : (" << mat.m() << ", " << mat.n() << ')' << nl
            << "Matrix columns != Vector size (" << x.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    auto tresult = tmp<Field<Type>>::New(mat.m(), Zero);
    auto& result = tresult.ref();

    for (label i = 0; i < mat.m(); ++i)
    {
        for (label j = 0; j < mat.n(); ++j)
        {
            result[i] += mat(i, j)*x[j];
        }
    }

    return tresult;
}


template<class Form, class Type>
template<class ListType>
Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::TmulImpl
(
    const ListType& x
) const
{
    const Matrix<Form, Type>& mat = *this;

    #ifdef FULLDEBUG
    if (mat.m() != x.size())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible Matrix and Vector:" << nl
            << "Matrix : (" << mat.m() << ", " << mat.n() << ')' << nl
            << "Matrix rows != Vector size (" << x.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    auto tresult = tmp<Field<Type>>::New(mat.n(), Zero);
    auto& result = tresult.ref();

    for (label i = 0; i < mat.m(); ++i)
    {
        const Type& val = x[i];
        for (label j = 0; j < mat.n(); ++j)
        {
            result[j] += val*mat(i, j);
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    checkSize();

    doAlloc();
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const Foam::zero)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    checkSize();

    doAlloc();

    std::fill_n(begin(), size(), Zero);
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const Type& val)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    checkSize();

    doAlloc();

    std::fill_n(begin(), size(), val);
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form, Type>& mat)
:
    mRows_(mat.mRows_),
    nCols_(mat.nCols_),
    v_(nullptr)
{
    if (mat.cdata())
    {
        doAlloc();

        std::copy(mat.cbegin(), mat.cend(), v_);
    }
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(Matrix<Form, Type>&& mat)
:
    mRows_(mat.mRows_),
    nCols_(mat.nCols_),
    v_(mat.v_)
{
    mat.mRows_ = 0;
    mat.nCols_ = 0;
    mat.v_ = nullptr;
}


template<class Form, class Type>
template<class Form2>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form2, Type>& mat)
:
    mRows_(mat.m()),
    nCols_(mat.n()),
    v_(nullptr)
{
    if (mat.cdata())
    {
        doAlloc();

        std::copy(mat.cbegin(), mat.cend(), v_);
    }
}


template<class Form, class Type>
template<class MatrixType>
inline Foam::Matrix<Form, Type>::Matrix
(
    const ConstMatrixBlock<MatrixType>& Mb
)
:
    mRows_(Mb.m()),
    nCols_(Mb.n()),
    v_(nullptr)
{
    doAlloc();

    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i,j);
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
inline Foam::Matrix<Form, Type>::Matrix
(
    const MatrixBlock<MatrixType>& Mb
)
:
    mRows_(Mb.m()),
    nCols_(Mb.n()),
    v_(nullptr)
{
    doAlloc();

    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i, j);
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::~Matrix()
{
    delete[] v_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::clear()
{
    if (v_)
    {
        delete[] v_;
        v_ = nullptr;
    }

    mRows_ = 0;
    nCols_ = 0;
}


template<class Form, class Type>
Foam::List<Type> Foam::Matrix<Form, Type>::release()
{
    List<Type> list;

    const label len = size();

    if (v_ && len)
    {
        UList<Type> storage(v_, len);
        list.swap(storage);

        v_ = nullptr;
    }
    clear();

    return list;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::swap(Matrix<Form, Type>& mat)
{
    if (this == &mat)
    {
        return;  // Self-swap is a no-op
    }

    std::swap(mRows_, mat.mRows_);
    std::swap(nCols_, mat.nCols_);
    std::swap(v_, mat.v_);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::transfer(Matrix<Form, Type>& mat)
{
    if (this == &mat)
    {
        return;  // Self-assignment is a no-op
    }

    clear();

    mRows_ = mat.mRows_;
    nCols_ = mat.nCols_;
    v_ = mat.v_;

    mat.mRows_ = 0;
    mat.nCols_ = 0;
    mat.v_ = nullptr;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::resize(const label m, const label n)
{
    if (m == mRows_ && n == nCols_)
    {
        return;
    }

    Matrix<Form, Type> newMatrix(m, n, Zero);

    const label mrow = min(m, mRows_);
    const label ncol = min(n, nCols_);

    for (label i = 0; i < mrow; ++i)
    {
        for (label j = 0; j < ncol; ++j)
        {
            newMatrix(i, j) = (*this)(i, j);
        }
    }

    transfer(newMatrix);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::resize_nocopy(const label mrow, const label ncol)
{
    if (mrow == mRows_ && ncol == nCols_)
    {
        return;
    }

    const label oldLen = (mRows_ * nCols_);

    const label newLen = (mrow * ncol);

    if (oldLen == newLen)
    {
        // Shallow resize is enough
        mRows_ = mrow;
        nCols_ = ncol;
    }
    else
    {
        this->clear();

        mRows_ = mrow;
        nCols_ = ncol;

        this->doAlloc();
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::round(const scalar tol)
{
    for (Type& val : *this)
    {
        if (mag(val) < tol)
        {
            val = Zero;
        }
    }
}


template<class Form, class Type>
Form Foam::Matrix<Form, Type>::T() const
{
    Form At(labelPair{n(), m()});

    for (label i = 0; i < m(); ++i)
    {
        for (label j = 0; j < n(); ++j)
        {
            At(j, i) = Detail::conj((*this)(i, j));
        }
    }

    return At;
}


template<class Form, class Type>
Form Foam::Matrix<Form, Type>::transpose() const
{
    Form At(labelPair{n(), m()});

    for (label i = 0; i < m(); ++i)
    {
        for (label j = 0; j < n(); ++j)
        {
            At(j, i) = (*this)(i, j);
        }
    }

    return At;
}


template<class Form, class Type>
Foam::List<Type> Foam::Matrix<Form, Type>::diag() const
{
    const label len = Foam::min(mRows_, nCols_);

    List<Type> result(len);

    for (label i=0; i < len; ++i)
    {
        result[i] = (*this)(i, i);
    }

    return result;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::diag(const UList<Type>& list)
{
    const label len = Foam::min(mRows_, nCols_);

    #ifdef FULLDEBUG
    if (list.size() != len)
    {
        FatalErrorInFunction
            << "List size (" << list.size()
            << ") incompatible with Matrix diagonal" << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        (*this)(i, i) = list[i];
    }
}


template<class Form, class Type>
Type Foam::Matrix<Form, Type>::trace() const
{
    const label len = Foam::min(mRows_, nCols_);

    Type val = Zero;

    for (label i=0; i < len; ++i)
    {
        val += (*this)(i, i);
    }

    return val;
}


template<class Form, class Type>
Foam::scalar Foam::Matrix<Form, Type>::columnNorm
(
    const label colIndex,
    const bool noSqrt
) const
{
    scalar result = Zero;

    for (label i=0; i < mRows_; ++i)
    {
        result += magSqr((*this)(i, colIndex));
    }

    return noSqrt ? result : Foam::sqrt(result);
}


template<class Form, class Type>
Foam::scalar Foam::Matrix<Form, Type>::norm(const bool noSqrt) const
{
    scalar result = Zero;

    for (const Type& val : *this)
    {
        result += magSqr(val);
    }

    return noSqrt ? result : Foam::sqrt(result);
}


template<class Form, class Type>
std::streamsize Foam::Matrix<Form, Type>::byteSize() const
{
    if (!is_contiguous<Type>::value)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << abort(FatalError);
    }
    return this->size_bytes();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Matrix<Form, Type>& mat)
{
    if (this == &mat)
    {
        return;  // Self-assignment is a no-op
    }

    if (mRows_ != mat.mRows_ || nCols_ != mat.nCols_)
    {
        clear();
        mRows_ = mat.mRows_;
        nCols_ = mat.nCols_;
        doAlloc();
    }

    if (v_)
    {
        std::copy(mat.cbegin(), mat.cend(), v_);
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(Matrix<Form, Type>&& mat)
{
    if (this != &mat)
    {
        // Self-assignment is a no-op
        this->transfer(mat);
    }
}


template<class Form, class Type>
template<class MatrixType>
void Foam::Matrix<Form, Type>::operator=
(
    const ConstMatrixBlock<MatrixType>& Mb
)
{
    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i, j);
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
void Foam::Matrix<Form, Type>::operator=
(
    const MatrixBlock<MatrixType>& Mb
)
{
    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i, j);
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Type& val)
{
    std::fill_n(begin(), size(), val);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Foam::zero)
{
    std::fill_n(begin(), size(), Zero);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator+=(const Matrix<Form, Type>& other)
{
    #ifdef FULLDEBUG
    if (this == &other)
    {
        FatalErrorInFunction
            << "Attempted addition to self"
            << abort(FatalError);
    }

    if (m() != other.m() || n() != other.n())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different sizes: ("
            << m() << ", " << n() << ") != ("
            << other.m() << ", " << other.n() << ')' << nl
            << abort(FatalError);
    }
    #endif

    auto iter2 = other.cbegin();
    for (Type& val : *this)
    {
        val += *iter2;
        ++iter2;
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator-=(const Matrix<Form, Type>& other)
{
    #ifdef FULLDEBUG
    if (this == &other)
    {
        FatalErrorInFunction
            << "Attempted subtraction from self"
            << abort(FatalError);
    }

    if (m() != other.m() || n() != other.n())
    {
        FatalErrorInFunction
            << "Attempt to subtract matrices with different sizes: ("
            << m() << ", " << n() << ") != ("
            << other.m() << ", " << other.n() << ')' << nl
            << abort(FatalError);
    }
    #endif

    auto iter2 = other.cbegin();
    for (Type& val : *this)
    {
        val -= *iter2;
        ++iter2;
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator+=(const Type& s)
{
    for (Type& val : *this)
    {
        val += s;
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator-=(const Type& s)
{
    for (Type& val : *this)
    {
        val -= s;
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator*=(const Type& s)
{
    for (Type& val : *this)
    {
        val *= s;
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator/=(const Type& s)
{
    for (Type& val : *this)
    {
        val /= s;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions * * * * * * * * * * * * * * //

//- Find max value in Matrix
template<class Form, class Type>
const Type& max(const Matrix<Form, Type>& mat)
{
    if (mat.empty())
    {
        FatalErrorInFunction
            << "Matrix is empty" << abort(FatalError);
    }

    return *(std::max_element(mat.cbegin(), mat.cend()));
}


//- Find min value in Matrix
template<class Form, class Type>
const Type& min(const Matrix<Form, Type>& mat)
{
    if (mat.empty())
    {
        FatalErrorInFunction
            << "Matrix is empty" << abort(FatalError);
    }

    return *(std::min_element(mat.cbegin(), mat.cend()));
}


//- Find the min/max values of Matrix
template<class Form, class Type>
MinMax<Type> minMax(const Matrix<Form, Type>& mat)
{
    MinMax<Type> result;

    for (const Type& val : mat)
    {
        result += val;
    }

    return result;
}



// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Matrix negation
template<class Form, class Type>
Form operator-(const Matrix<Form, Type>& mat)
{
    Form result(mat.sizes());

    std::transform
    (
        mat.cbegin(),
        mat.cend(),
        result.begin(),
        std::negate<Type>()
    );

    return result;
}


//- Matrix addition. Returns Matrix of the same form as the first parameter.
template<class Form1, class Form2, class Type>
Form1 operator+
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B
)
{
    #ifdef FULLDEBUG
    if (A.m() != B.m() || A.n() != B.n())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different sizes: ("
            << A.m() << ", " << A.n() << ") != ("
            << B.m() << ", " << B.n() << ')' << nl
            << abort(FatalError);
    }
    #endif

    Form1 result(A.sizes());

    std::transform
    (
        A.cbegin(),
        A.cend(),
        B.cbegin(),
        result.begin(),
        std::plus<Type>()
    );

    return result;
}


//- Matrix subtraction. Returns Matrix of the same form as the first parameter.
template<class Form1, class Form2, class Type>
Form1 operator-
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B
)
{
    #ifdef FULLDEBUG
    if (A.m() != B.m() || A.n() != B.n())
    {
        FatalErrorInFunction
            << "Attempt to subtract matrices with different sizes: ("
            << A.m() << ", " << A.n() << ") != ("
            << B.m() << ", " << B.n() << ')' << nl
            << abort(FatalError);
    }
    #endif

    Form1 result(A.sizes());

    std::transform
    (
        A.cbegin(),
        A.cend(),
        B.cbegin(),
        result.begin(),
        std::minus<Type>()
    );

    return result;
}


//- Scalar multiplication of Matrix
template<class Form, class Type>
Form operator*(const Type& s, const Matrix<Form, Type>& mat)
{
    Form result(mat.sizes());

    std::transform
    (
        mat.cbegin(),
        mat.cend(),
        result.begin(),
        [&](const Type& val) { return s * val; }
    );

    return result;
}


//- Scalar multiplication of Matrix
template<class Form, class Type>
Form operator*(const Matrix<Form, Type>& mat, const Type& s)
{
    return s*mat;
}


//- Scalar addition of Matrix
template<class Form, class Type>
Form operator+(const Type& s, const Matrix<Form, Type>& mat)
{
    Form result(mat.sizes());

    std::transform
    (
        mat.cbegin(),
        mat.cend(),
        result.begin(),
        [&](const Type& val) { return s + val; }
    );

    return result;
}


//- Scalar addition of Matrix
template<class Form, class Type>
Form operator+(const Matrix<Form, Type>& mat, const Type& s)
{
    return s + mat;
}


//- Scalar subtraction of Matrix
template<class Form, class Type>
Form operator-(const Type& s, const Matrix<Form, Type>& mat)
{
    Form result(mat.sizes());

    std::transform
    (
        mat.cbegin(),
        mat.end(),
        result.begin(),
        [&](const Type& val) { return s - val; }
    );

    return result;
}


//- Scalar subtraction of Matrix
template<class Form, class Type>
Form operator-(const Matrix<Form, Type>& mat, const Type& s)
{
    Form result(mat.sizes());

    std::transform
    (
        mat.cbegin(),
        mat.end(),
        result.begin(),
        [&](const Type& val) { return val - s; }
    );

    return result;
}


//- Scalar division of Matrix
template<class Form, class Type>
Form operator/(const Matrix<Form, Type>& mat, const Type& s)
{
    Form result(mat.sizes());

    std::transform
    (
        mat.cbegin(),
        mat.end(),
        result.begin(),
        [&](const Type& val) { return val / s; }
    );

    return result;
}


//- Matrix-Matrix multiplication using ikj-algorithm
template<class Form1, class Form2, class Type>
typename typeOfInnerProduct<Type, Form1, Form2>::type
operator*
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B
)
{
    #ifdef FULLDEBUG
    if (A.n() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible matrices:" << nl
            << "Matrix A : (" << A.m() << ", " << A.n() << ')' << nl
            << "Matrix B : (" << B.m() << ", " << B.n() << ')' << nl
            << "The columns of A must equal rows of B"
            << abort(FatalError);
    }
    #endif

    typename typeOfInnerProduct<Type, Form1, Form2>::type AB
    (
        A.m(),
        B.n(),
        Zero
    );

    for (label i = 0; i < AB.m(); ++i)
    {
        for (label k = 0; k < B.m(); ++k)
        {
            for (label j = 0; j < AB.n(); ++j)
            {
                AB(i, j) += A(i, k)*B(k, j);
            }
        }
    }

    return AB;
}


//- Implicit inner product of Matrix-Matrix, equivalent to A.T()*B
template<class Form1, class Form2, class Type>
typename typeOfInnerProduct<Type, Form1, Form2>::type
operator&
(
    const Matrix<Form1, Type>& AT,
    const Matrix<Form2, Type>& B
)
{
    #ifdef FULLDEBUG
    if (AT.m() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible matrices:" << nl
            << "Matrix A : (" << AT.m() << ", " << AT.n() << ')' << nl
            << "Matrix B : (" << B.m() << ", " << B.n() << ')' << nl
            << "The rows of A must equal rows of B"
            << abort(FatalError);
    }
    #endif

    typename typeOfInnerProduct<Type, Form1, Form2>::type AB
    (
        AT.n(),
        B.n(),
        Zero
    );

    for (label k = 0; k < B.m(); ++k)
    {
        for (label i = 0; i < AB.m(); ++i)
        {
            for (label j = 0; j < AB.n(); ++j)
            {
                AB(i, j) += Detail::conj(AT(k, i))*B(k, j);
            }
        }
    }

    return AB;
}


//- Implicit outer product of Matrix-Matrix, equivalent to A*B.T()
template<class Form1, class Form2, class Type>
typename typeOfInnerProduct<Type, Form1, Form2>::type
operator^
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& BT
)
{
    #ifdef FULLDEBUG
    if (A.n() != BT.n())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible matrices:" << nl
            << "Matrix A : (" << A.m() << ", " << A.n() << ')' << nl
            << "Matrix B : (" << BT.m() << ", " << BT.n() << ')' << nl
            << "The columns of A must equal columns of B"
            << abort(FatalError);
    }
    #endif

    typename typeOfInnerProduct<Type, Form1, Form2>::type AB
    (
        A.m(),
        BT.m(),
        Zero
    );

    for (label i = 0; i < AB.m(); ++i)
    {
        for (label j = 0; j < AB.n(); ++j)
        {
            for (label k = 0; k < BT.n(); ++k)
            {
                AB(i, j) += A(i, k)*Detail::conj(BT(j, k));
            }
        }
    }

    return AB;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
