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

#include "MatrixBlock.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MatrixType>
Foam::ConstMatrixBlock<MatrixType>::operator Field<cmptType>() const
{
    if (nCols_ != 1)
    {
        FatalErrorInFunction
            << "Number of columns " << nCols_ << " != 1"
            << abort(FatalError);
    }

    Field<cmptType> f(mRows_);

    forAll(f, i)
    {
        f[i] = operator()(i, 0);
    }

    return f;
}


template<class MatrixType>
Foam::MatrixBlock<MatrixType>::operator Field<cmptType>() const
{
    if (nCols_ != 1)
    {
        FatalErrorInFunction
            << "Number of columns " << nCols_ << " != 1"
            << abort(FatalError);
    }

    Field<cmptType> f(mRows_);

    forAll(f, i)
    {
        f[i] = operator()(i, 0);
    }

    return f;
}


template<class MatrixType>
Foam::label Foam::ConstMatrixBlock<MatrixType>::disallow
(
    const char* what
) const
{
    FatalErrorInFunction
        << "Block addresses " << what
        << " outside matrix or invalid matrix components"
        << abort(FatalError);
    return 0;
}


template<class MatrixType>
Foam::label Foam::MatrixBlock<MatrixType>::disallow
(
    const char* what
) const
{
    FatalErrorInFunction
        << "Block addresses " << what
        << " outside matrix or invalid matrix components"
        << abort(FatalError);
    return 0;
}


template<class MatrixType> void Foam::ConstMatrixBlock<MatrixType>::checkIndex
(
    const label i,
    const label j
) const
{
    if (i < 0 || i >= mRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " is out of range 0 ... " << mRows_ - 1
            << abort(FatalError);
    }
    else if (j < 0 || j >= nCols_)
    {
        FatalErrorInFunction
            << "Index " << j << " is out of range 0 ... " << nCols_ - 1
            << abort(FatalError);
    }
}


template<class MatrixType> void Foam::MatrixBlock<MatrixType>::checkIndex
(
    const label i,
    const label j
) const
{
    if (i < 0 || i >= mRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " is out of range 0 ... " << mRows_ - 1
            << abort(FatalError);
    }
    else if (j < 0 || j >= nCols_)
    {
        FatalErrorInFunction
            << "Index " << j << " is out of range 0 ... " << nCols_ - 1
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class MatrixType>
template<class Form>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const Matrix<Form, cmptType>& Mb
)
{
    if (mRows_ != Mb.m() || nCols_ != Mb.n())
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << mRows_ << "x" << nCols_ << " != "
            << Mb.m() << "x" << Mb.n()
            << abort(FatalError);
    }

    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i, j);
        }
    }
}


template<class MatrixType>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const ConstMatrixBlock<MatrixType>& Mb
)
{
    if (reinterpret_cast<const ConstMatrixBlock<MatrixType>*>(this) != &Mb)
    {
        if (mRows_ != Mb.m() || nCols_ != Mb.n())
        {
            FatalErrorInFunction
                << "Attempt to assign blocks of different sizes: "
                << mRows_ << "x" << nCols_ << " != "
                << Mb.m() << "x" << Mb.n()
                << abort(FatalError);
        }

        for (label i = 0; i < mRows_; ++i)
        {
            for (label j = 0; j < nCols_; ++j)
            {
                (*this)(i, j) = Mb(i, j);
            }
        }
    }
}


template<class MatrixType>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const MatrixBlock<MatrixType>& Mb
)
{
    if (this != &Mb)
    {
        if (mRows_ != Mb.m() || nCols_ != Mb.n())
        {
            FatalErrorInFunction
                << "Attempt to assign blocks of different sizes: "
                << mRows_ << "x" << nCols_ << " != "
                << Mb.m() << "x" << Mb.n()
                << abort(FatalError);
        }

        for (label i = 0; i < mRows_; ++i)
        {
            for (label j = 0; j < nCols_; ++j)
            {
                (*this)(i, j) = Mb(i, j);
            }
        }
    }
}


template<class MatrixType>
template<class MatrixType2>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const ConstMatrixBlock<MatrixType2>& Mb
)
{
    if (reinterpret_cast<const ConstMatrixBlock<MatrixType2>*>(this) != &Mb)
    {
        if (mRows_ != Mb.m() || nCols_ != Mb.n())
        {
            FatalErrorInFunction
                << "Attempt to assign blocks of different sizes: "
                << mRows_ << "x" << nCols_ << " != "
                << Mb.m() << "x" << Mb.n()
                << abort(FatalError);
        }

        for (label i = 0; i < mRows_; ++i)
        {
            for (label j = 0; j < nCols_; ++j)
            {
                (*this)(i, j) = Mb(i, j);
            }
        }
    }
}


template<class MatrixType>
template<class MatrixType2>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const MatrixBlock<MatrixType2>& Mb
)
{
    if (this != &Mb)
    {
        if (mRows_ != Mb.m() || nCols_ != Mb.n())
        {
            FatalErrorInFunction
                << "Attempt to assign blocks of different sizes: "
                << mRows_ << "x" << nCols_ << " != "
                << Mb.m() << "x" << Mb.n()
                << abort(FatalError);
        }

        for (label i = 0; i < mRows_; ++i)
        {
            for (label j = 0; j < nCols_; ++j)
            {
                (*this)(i, j) = Mb(i, j);
            }
        }
    }
}


template<class MatrixType>
template
<
    template<class, Foam::direction, Foam::direction> class MSBlock,
    class SubTensor,
    Foam::direction BRowStart,
    Foam::direction BColStart
>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const MSBlock<SubTensor, BRowStart, BColStart>& Mb
)
{
    if (mRows_ != Mb.mRows || nCols_ != Mb.nCols)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << mRows_ << "x" << nCols_ << " != "
            << Mb.mRows << "x" << Mb.nCols
            << abort(FatalError);
    }

    for (direction i = 0; i < mRows_; ++i)
    {
        for (direction j = 0; j < nCols_; ++j)
        {
            operator()(i, j) = Mb(i, j);
        }
    }
}


template<class MatrixType>
template
<
    template<class, Foam::direction> class VSBlock,
    class SubVector,
    Foam::direction BStart
>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const VSBlock<SubVector, BStart>& Mb
)
{
    if (mRows_ != Mb.nComponents || nCols_ != 1)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << mRows_ << "x" << nCols_ << " != "
            << Mb.nComponents << "x" << 1
            << abort(FatalError);
    }

    for (direction i = 0; i < mRows_; ++i)
    {
        operator()(i, 0) = Mb[i];
    }
}


template<class MatrixType>
template<class MSForm, Foam::direction Nrows, Foam::direction Ncols>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const MatrixSpace<MSForm, cmptType, Nrows, Ncols>& ms
)
{
    if (mRows_ != Nrows || nCols_ != Ncols)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << mRows_ << "x" << nCols_ << " != "
            << Nrows << "x" << Ncols
            << abort(FatalError);
    }

    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = ms(i, j);
        }
    }
}


template<class MatrixType>
template<class VSForm, Foam::direction Ncmpts>
void Foam::MatrixBlock<MatrixType>::operator=
(
    const VectorSpace<VSForm, cmptType, Ncmpts>& ms
)
{
    if (mRows_ != Ncmpts || nCols_ != 1)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << mRows_ << "x" << nCols_ << " != "
            << Ncmpts << "x" << 1
            << abort(FatalError);
    }

    for (direction i = 0; i < Ncmpts; ++i)
    {
        operator()(i, 0) = ms[i];
    }
}


template<class MatrixType>
void Foam::MatrixBlock<MatrixType>::operator=(const Field<cmptType>& f)
{
    if (mRows_ != f.size() || nCols_ != 1)
    {
        FatalErrorInFunction
            << "Error: cannot assign blocks of different size (left is "
            << mRows_ << "x" << nCols_ << " != "
            << f.size() << "x" << 1
            << abort(FatalError);
    }

    forAll(f, i)
    {
        operator()(i, 0) = f[i];
    }
}


// ************************************************************************* //
