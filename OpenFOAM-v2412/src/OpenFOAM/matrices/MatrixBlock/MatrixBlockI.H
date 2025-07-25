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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MatrixType>
Foam::ConstMatrixBlock<MatrixType>::ConstMatrixBlock
(
    const MatrixType& matrix,
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
:
    matrix_(matrix),
    mRows_(0 < m ? m : disallow("row dim")),
    nCols_(0 < n ? n : disallow("col dim")),
    rowStart_
    (
        0 <= mStart
     || mStart + mRows_ <= matrix.m()
     ?  mStart : disallow("row start")
    ),
    colStart_
    (
        0 <= nStart
     || nStart + nCols_ <= matrix.n()
     ?  nStart : disallow("col start")
    )
{}


template<class MatrixType>
Foam::MatrixBlock<MatrixType>::MatrixBlock
(
    MatrixType& matrix,
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
:
    matrix_(matrix),
    mRows_(0 < m ? m : disallow("row dim")),
    nCols_(0 < n ? n : disallow("col dim")),
    rowStart_
    (
        0 <= mStart
     || mStart + mRows_ <= matrix.m()
     ?  mStart : disallow("row start")
    ),
    colStart_
    (
        0 <= nStart
     || nStart + nCols_ <= matrix.n()
     ?  nStart : disallow("col start")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MatrixType>
inline Foam::label Foam::ConstMatrixBlock<MatrixType>::m() const
{
    return mRows_;
}


template<class MatrixType>
inline Foam::label Foam::ConstMatrixBlock<MatrixType>::n() const
{
    return nCols_;
}


template<class MatrixType>
inline Foam::label Foam::MatrixBlock<MatrixType>::m() const
{
    return mRows_;
}


template<class MatrixType>
inline Foam::label Foam::MatrixBlock<MatrixType>::n() const
{
    return nCols_;
}


template<class MatrixType>
inline Foam::labelPair Foam::ConstMatrixBlock<MatrixType>::sizes() const
{
    return labelPair(mRows_, nCols_);
}


template<class MatrixType>
inline Foam::labelPair Foam::MatrixBlock<MatrixType>::sizes() const
{
    return labelPair(mRows_, nCols_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class MatrixType>
inline const typename MatrixType::cmptType&
Foam::ConstMatrixBlock<MatrixType>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    checkIndex(i, j);
    #endif

    return matrix_(i + rowStart_, j + colStart_);
}


template<class MatrixType>
inline const typename MatrixType::cmptType&
Foam::MatrixBlock<MatrixType>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    checkIndex(i, j);
    #endif

    return matrix_(i + rowStart_, j + colStart_);
}


template<class MatrixType>
inline typename MatrixType::cmptType&
Foam::MatrixBlock<MatrixType>::operator()
(
    const label i,
    const label j
)
{
    #ifdef FULLDEBUG
    checkIndex(i, j);
    #endif

    return matrix_(i + rowStart_, j + colStart_);
}


// ************************************************************************* //
