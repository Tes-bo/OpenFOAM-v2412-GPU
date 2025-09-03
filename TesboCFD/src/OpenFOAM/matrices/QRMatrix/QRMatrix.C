/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "QRMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MatrixType>
Foam::label Foam::QRMatrix<MatrixType>::calcShapeFactor
(
    const MatrixType& A
) const
{
    if (mode_ == modes::FULL)
    {
        return A.m();
    }
    else if (mode_ == modes::ECONOMY)
    {
        return min(A.m(), A.n());
    }

    return 0;
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::decompose
(
    MatrixType& AT
)
{
    const label n = AT.m();
    const label m = AT.n();

    List<cmptType> Rdiag(n, Zero);

    for (label k = 0; k < n; ++k)
    {
        // Compute 2-norm of k-th column without under/overflow
        scalar nrm = 0;

        for (label i = k; i < m; ++i)
        {
            nrm = Foam::hypot(nrm, mag(AT(k,i)));
        }

        // Avoid divide-by-zero. Use compare with == 0 not mag or VSMALL etc,
        // otherwise wrong results (may need more investigation)
        if (nrm != scalar(0))
        {
            // Form k-th Householder vector
            if (mag(AT(k,k)) < 0)
            {
                nrm *= -1;
            }

            for (label i = k; i < m; ++i)
            {
                AT(k,i) /= nrm;
            }

            AT(k,k) += scalar(1);

            // Apply transformation to remaining columns
            for (label j = k + 1; j < n; ++j)
            {
                cmptType s(Zero);

                for (label i = k; i < m; ++i)
                {
                    s += Detail::conj(AT(k,i))*AT(j,i);
                }

                if (mag(AT(k,k)) > SMALL)
                {
                    s /= -AT(k,k);
                }

                for (label i = k; i < m; ++i)
                {
                    AT(j,i) += s*AT(k,i);
                }
            }
        }

        Rdiag[k] = -nrm;
    }

    calcQ(AT);

    calcR(AT, Rdiag);
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::decompose
(
    MatrixType& AT,
    const bool pivoting
)
{
    const label n = AT.m();
    const label m = AT.n();
    const label sz = min(m,n);

    // Initialise permutation vector, and column norm vector
    p_ = identity(n);

    // Initialise vector norms of each column of A
    List<scalar> norms(n, Zero);
    for (label k = 0; k < n; ++k)
    {
        for (label i = 0; i < m; ++i)
        {
            norms[k] += magSqr(AT(k, i));
        }
    }

    List<cmptType> Rdiag(n, Zero);

    for (label k = 0; k < sz; ++k)
    {
        const auto it = std::next(norms.cbegin(), k);
        const label maxNormi =
            std::distance
            (
                it,
                std::max_element(it, norms.cend())
            );

        // Swap A, P and norms according to pivot column if the current
        // column is not the max norm column. Also avoid any column swaps
        // where the leading element is very small
        if (mag(AT(k + maxNormi,k)) > SMALL && maxNormi != 0)
        {
            const RMatrix R1(AT.subRow(k));
            const RMatrix R2(AT.subRow(maxNormi + k));

            AT.subRow(k) = R2;
            AT.subRow(maxNormi + k) = R1;

            std::swap(p_[k], p_[maxNormi + k]);
            std::swap(norms[k], norms[maxNormi + k]);
        }

        {
            // Compute 2-norm of k-th column without under/overflow
            scalar nrm = 0;

            for (label i = k; i < m; ++i)
            {
                nrm = Foam::hypot(nrm, mag(AT(k,i)));
            }

            if (nrm != scalar(0))
            {
                // Form k-th Householder vector
                if (mag(AT(k,k)) < 0)
                {
                    nrm *= -1;
                }

                for (label i = k; i < m; ++i)
                {
                    AT(k,i) /= nrm;
                }

                AT(k,k) += scalar(1);

                // Apply transformation to remaining columns
                for (label j = k + 1; j < n; ++j)
                {
                    cmptType s(Zero);

                    for (label i = k; i < m; ++i)
                    {
                        s += Detail::conj(AT(k,i))*AT(j,i);
                    }

                    if (mag(AT(k,k)) > SMALL)
                    {
                        s /= -AT(k,k);
                    }

                    for (label i = k; i < m; ++i)
                    {
                        AT(j,i) += s*AT(k,i);
                    }
                }
            }

            Rdiag[k] = -nrm;
        }

        // Update norms
        if (k < sz - 1)
        {
            label q = k + 1;
            for (const auto& val : RMatrix(AT.subColumn(k, q)))
            {
                norms[q] -= magSqr(val);
                ++q;
            }
        }
    }

    calcQ(AT);

    calcR(AT, Rdiag);
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::calcQ
(
    const MatrixType& AT
)
{
    if (output_ == outputs::ONLY_R)
    {
        return;
    }

    const label n = AT.m();
    const label m = AT.n();

    Q_.resize(m, sz_);

    MatrixType QT(Q_.transpose());

    for (label k = sz_ - 1; k >= 0; --k)
    {
        for (label i = 0; i < m; ++i)
        {
            QT(k,i) = Zero;
        }

        QT(k,k) = scalar(1);

        for (label j = k; j < sz_; ++j)
        {
            if (n > k && mag(AT(k,k)) > SMALL)
            {
                cmptType s(Zero);

                for (label i = k; i < m; ++i)
                {
                    s += AT(k,i)*QT(j,i);
                }

                s /= -AT(k,k);

                for (label i = k; i < m; ++i)
                {
                    QT(j,i) += s*Detail::conj(AT(k,i));
                }
            }
        }
    }

    Q_ = QT.T();
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::calcR
(
    const MatrixType& AT,
    const List<cmptType>& diag
)
{
    if (output_ == outputs::ONLY_Q)
    {
        return;
    }

    const label n = AT.m();

    R_.resize(sz_, n);

    for (label i = 0; i < sz_; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            if (i < j)
            {
                R_(i,j) = AT(j,i);
            }
            else if (i == j)
            {
                R_(i,j) = diag[i];
            }
        }
    }
}


template<class MatrixType>
template<template<typename> class ListContainer>
void Foam::QRMatrix<MatrixType>::solvex
(
    ListContainer<cmptType>& x
) const
{
    const label n = R_.n();

    for (label i = n - 1; 0 <= i; --i)
    {
        cmptType sum = x[i];

        for (label j = i + 1; j < n; ++j)
        {
            sum -= x[j]*R_(i, j);
        }

        if (mag(R_(i, i)) < SMALL)
        {
            FatalErrorInFunction
                << "Back-substitution failed due to small diagonal"
                << abort(FatalError);
        }

        x[i] = sum/R_(i, i);
    }
}


template<class MatrixType>
template<template<typename> class ListContainer>
void Foam::QRMatrix<MatrixType>::solveImpl
(
    List<cmptType>& x,
    const ListContainer<cmptType>& source
) const
{
    // Assert (&x != &source) ?
    const label m = Q_.m();

    // x = Q_.T()*source;  i.e., Q_.Tmul(source)
    for (label i = 0; i < m; ++i)
    {
        x[i] = 0;
        for (label j = 0; j < m; ++j)
        {
            x[i] += Q_(j, i)*source[j];
        }
    }

    solvex(x);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MatrixType>
Foam::QRMatrix<MatrixType>::QRMatrix
(
    const modes mode,
    const outputs output,
    const bool pivoting,
    MatrixType& A
)
:
    mode_(mode),
    output_(output),
    sz_(calcShapeFactor(A)),
    Q_(),
    R_(),
    p_()
{
    // Non-conjugate transpose of input matrix
    MatrixType AT(A.transpose());
    A.clear();

    if (pivoting)
    {
        decompose(AT, pivoting);
    }
    else
    {
        decompose(AT);
    }

    AT.clear();
}


template<class MatrixType>
Foam::QRMatrix<MatrixType>::QRMatrix
(
    const MatrixType& A,
    const modes mode,
    const outputs output,
    const bool pivoting
)
:
    mode_(mode),
    output_(output),
    sz_(calcShapeFactor(A)),
    Q_(),
    R_(),
    p_()
{
    MatrixType AT(A.transpose());

    if (pivoting)
    {
        decompose(AT, pivoting);
    }
    else
    {
        decompose(AT);
    }

    AT.clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class MatrixType>
void Foam::QRMatrix<MatrixType>::solve
(
    List<cmptType>& x,
    const UList<cmptType>& source
) const
{
    solveImpl(x, source);
}


template<class MatrixType>
template<class Addr>
void Foam::QRMatrix<MatrixType>::solve
(
    List<cmptType>& x,
    const IndirectListBase<cmptType, Addr>& source
) const
{
    solveImpl(x, source);
}


template<class MatrixType>
Foam::tmp<Foam::Field<typename MatrixType::cmptType>>
Foam::QRMatrix<MatrixType>::solve
(
    const UList<cmptType>& source
) const
{
    auto tresult(Q_.Tmul(source));

    solvex(tresult.ref());

    return tresult;
}


template<class MatrixType>
template<class Addr>
Foam::tmp<Foam::Field<typename MatrixType::cmptType>>
Foam::QRMatrix<MatrixType>::solve
(
    const IndirectListBase<cmptType, Addr>& source
) const
{
    auto tresult(Q_.Tmul(source));

    solvex(tresult.ref());

    return tresult;
}


template<class MatrixType>
Foam::RectangularMatrix<typename MatrixType::cmptType>
Foam::QRMatrix<MatrixType>::solve
(
    const RMatrix& rhs
)
{
    const label mRows = R_.m();
    const label nCols = R_.n();
    const label pCols = rhs.n();

    #ifdef FULLDEBUG
    {
        const label qRows = rhs.m();
        if (mRows != qRows)
        {
            FatalErrorInFunction
                << "Linear system is not solvable since the number of rows of"
                << "A and rhs are not equal:" << tab << mRows << "vs." << qRows
                << abort(FatalError);
        }

        const List<cmptType> diag(R_.diag());
        for (const auto& val : diag)
        {
            if (mag(val) < SMALL)
            {
                WarningInFunction
                    << "SMALL-valued diagonal elem in back-substitution."
                    << endl;
            }
        }
    }
    #endif

    RMatrix b({nCols, pCols}, Zero);

    for (label i = 0; i < pCols; ++i)
    {
        for (label j = mRows - 1; -1 < j; --j)
        {
            cmptType alpha(rhs(j, i));

            for (label k = j + 1; k < mRows; ++k)
            {
                alpha -= b(k, i)*R_(j, k);
            }

            b(j, i) = alpha/R_(j, j);
        }
    }

    return b;
}


template<class MatrixType>
typename Foam::QRMatrix<MatrixType>::SMatrix
Foam::QRMatrix<MatrixType>::inv() const
{
    const label m = Q_.m();

    Field<cmptType> x(m);
    SMatrix inv(m);

    for (label i = 0; i < m; ++i)
    {
        for (label j = 0; j < m; ++j)
        {
            x[j] = Q_(i, j);
        }
        solvex(x);
        inv.subColumn(i) = x;
    }

    return inv;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class MatrixType>
MatrixType Foam::MatrixTools::pinv
(
    const MatrixType& A,
    scalar tol
)
{
    typedef typename MatrixType::cmptType cmptType;

    if (A.empty())
    {
        FatalErrorInFunction
            << "Empty matrix found."
            << abort(FatalError);
    }

    if (A.size() == 1)
    {
        if (A(0,0) == cmptType(0))
        {
            return MatrixType({1, 1}, cmptType(0));
        }

        return MatrixType({1, 1}, cmptType(1)/(A(0,0) + cmptType(VSMALL)));
    }

    QRMatrix<MatrixType> QRM
    (
        A,
        QRMatrix<MatrixType>::modes::FULL,
        QRMatrix<MatrixType>::outputs::BOTH_QR,
        QRMatrix<MatrixType>::pivoting::TRUE
    );
    const MatrixType& R = QRM.R();
    const MatrixType& Q = QRM.Q();

    // R1 (KP:p. 648)
    // Find the first diagonal element index with (almost) zero value in R
    label elemi = 0;
    label i = 0;
    while (i < 2)
    {
        const List<cmptType> diag(R.diag());

        auto lessThan = [=](const cmptType& x) { return tol > mag(x); };

        elemi =
            std::distance
            (
                diag.cbegin(),
                std::find_if(diag.cbegin(), diag.cend(), lessThan)
            );

        if (elemi == 0)
        {
            if (i == 0)
            {
                WarningInFunction
                    << "The largest diagonal element magnitude is nearly zero. "
                    << "Tightening the tolerance."
                    << endl;

                tol = 1e-13;
                ++i;
            }
            else
            {
                WarningInFunction
                    << "The largest diagonal element magnitude is nearly zero. "
                    << "Returning a zero matrix."
                    << endl;
                ++i;

                return MatrixType({A.n(), A.m()}, Zero);
            }
        }
        else
        {
            i += 2;
        }
    }

    // Exclude the first (almost) zero diagonal row and the rows below
    // since R diagonal is already descending due to the QR column pivoting
    const RectangularMatrix<cmptType> R1(R.subMatrix(0, 0, elemi));

    // R2 (KP:p. 648)
    if (R1.n() < R1.m())
    {
        const SquareMatrix<cmptType> C(R1&R1);

        QRMatrix<SquareMatrix<cmptType>> QRSolve
        (
            C,
            QRMatrix<SquareMatrix<cmptType>>::modes::FULL,
            QRMatrix<SquareMatrix<cmptType>>::outputs::BOTH_QR
        );

        RectangularMatrix<cmptType> R2
        (
            QRSolve.solve
            (
                QRSolve.Q() & R1.T()
            )
        );

        // R3 (KP:p. 648)
        R2.resize(R.m(), R.n());

        return MatrixType((QRM.P()^R2)^Q);
    }
    else
    {
        const SquareMatrix<cmptType> C(R1^R1);

        QRMatrix<SquareMatrix<cmptType>> QRSolve
        (
            C,
            QRMatrix<SquareMatrix<cmptType>>::modes::FULL,
            QRMatrix<SquareMatrix<cmptType>>::outputs::BOTH_QR
        );

        RectangularMatrix<cmptType> R2
        (
            QRSolve.solve
            (
                QRSolve.Q() & R1
            )
        );

        // R3
        R2.resize(R.m(), R.n());

        return MatrixType((QRM.P()^R2)^Q);
    }
}


// ************************************************************************* //
