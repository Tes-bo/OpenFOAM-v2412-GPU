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

Class
    Foam::QRMatrix

Description
    \c QRMatrix computes QR decomposition of a given
    scalar/complex matrix \c A into the following:

    \verbatim
        A = Q R
    \endverbatim

    or in case of QR decomposition with column pivoting:

    \verbatim
        A P = Q R
    \endverbatim

    where
    \vartable
        Q | Unitary/orthogonal matrix
        R | Upper triangular matrix
        P | Permutation matrix
    \endvartable

    References:
    \verbatim
        TNT implementation:
            Pozo, R. (1997).
            Template Numerical Toolkit for linear algebra:
            High performance programming with C++
            and the Standard Template Library.
            The International Journal of Supercomputer Applications
            and High Performance Computing, 11(3), 251-263.
            DOI:10.1177/109434209701100307

        QR decomposition with column pivoting (tag:QSB):
            Quintana-Ortí, G., Sun, X., & Bischof, C. H. (1998).
            A BLAS-3 version of the QR factorization with column pivoting.
            SIAM Journal on Scientific Computing, 19(5), 1486-1494.
            DOI:10.1137/S1064827595296732

        Moore-Penrose inverse algorithm (tags:KP; KPP):
            Katsikis, V. N., & Pappas, D. (2008).
            Fast computing of the Moore-Penrose inverse matrix.
            Electronic Journal of Linear Algebra, 17(1), 637-650.
            DOI:10.13001/1081-3810.1287

            Katsikis, V. N., Pappas, D., & Petralias, A. (2011).
            An improved method for the computation of
            the Moore–Penrose inverse matrix.
            Applied Mathematics and Computation, 217(23), 9828-9834.
            DOI:10.1016/j.amc.2011.04.080

        Tolerance for the Moore-Penrose inverse algorithm (tag:TA):
            Toutounian, F., & Ataei, A. (2009).
            A new method for computing Moore–Penrose inverse matrices.
            Journal of Computational and applied Mathematics, 228(1), 412-417.
            DOI:10.1016/j.cam.2008.10.008
    \endverbatim

Usage

    \heading Input:
    \vartable
        A | RectangularMatrix<Type> or SquareMatrix<Type>
    \endvartable

    \heading Options for the decomposition mode:
    \vartable
        modes::FULL     | compute full-size decomposition
        modes::ECONOMY  | compute economy-size decomposition
    \endvartable

    \heading Options for the output types:
    \vartable
        outputs::ONLY\_Q     | compute only Q
        outputs::ONLY\_R     | compute only R
        outputs::BOTH\_QR    | compute both Q and R
    \endvartable

    \heading Options for the column pivoting:
    \vartable
        pivoting::FALSE     | switch off column pivoting
        pivoting::TRUE      | switch on column pivoting
    \endvartable

    \heading Output:
    \vartable
        Q | m-by-m (FULL) or m-by-k (ECONOMY) with k = min(m,n)
        R | m-by-n (FULL) or k-by-n (ECONOMY) with k = min(m,n)
        p | n-element label list
        P | n-by-n permutation matrix
    \endvartable

Notes
    - \c QRMatrix involves modified implementations of the public-domain
    library \c TNT, which is available at https://math.nist.gov/tnt/index.html.
    - \c Q and \c R are always of the same type of the matrix \c A.
    - \c Type can be \c scalar or \c complex.

See also
    Test-QRMatrix.C

SourceFiles
    QRMatrix.C
    QRMatrixI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_QRMatrix_H
#define Foam_QRMatrix_H

#include "RectangularMatrix.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class QRMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class QRMatrix
{
public:

    typedef typename MatrixType::cmptType cmptType;
    typedef SquareMatrix<cmptType> SMatrix;
    typedef RectangularMatrix<cmptType> RMatrix;

    //- Options for the decomposition mode
    enum modes : uint8_t
    {
        FULL = 1,         //!< compute full-size decomposition
        ECONOMY = 2,      //!< compute economy-size decomposition
    };

    //- Options for the output types
    enum outputs : uint8_t
    {
        ONLY_Q = 1,      //!< compute only Q
        ONLY_R = 2,      //!< compute only R
        BOTH_QR = 3      //!< compute both Q and R
    };

    //- Options for the column pivoting
    enum pivoting : bool
    {
        FALSE = false,      //!< switch off column pivoting
        TRUE  = true        //!< switch on column pivoting
    };


private:

    // Private Data

        //- Decomposition mode
        const modes mode_;

        //- Output type
        const outputs output_;

        //- Output shape factor based on the decomposition mode
        const label sz_;

        //- Unitary/orthogonal matrix
        MatrixType Q_;

        //- Upper triangular matrix
        MatrixType R_;

        //- Permutation list (if column-pivoting is on)
        labelList p_;


    // Private Member Functions

    // Evaluation

        //- Calculate output shape factor based on the decomposition mode
        label calcShapeFactor(const MatrixType& A) const;

        //- Compute QR decomposition
        void decompose(MatrixType& AT);

        //- Compute QR decomposition with column pivoting
        void decompose(MatrixType& AT, const bool pivot);

        //- Calculate Q
        void calcQ(const MatrixType& AT);

        //- Calculate R
        void calcR(const MatrixType& AT, const List<cmptType>& diag);


    // Linear system solution

        //- Solve the linear system with the Field argument x
        //- initialized to the appropriate transformed source
        //- (e.g. Q.T()*source) and return the solution in x
        template<template<typename> class ListContainer>
        void solvex(ListContainer<cmptType>& x) const;

        //- Solve the linear system with the given source
        //- and return the solution in x
        template<template<typename> class ListContainer>
        void solveImpl
        (
            List<cmptType>& x,
            const ListContainer<cmptType>& source
        ) const;


        //- No copy construct
        QRMatrix(const QRMatrix&) = delete;

        //- No copy assignment
        QRMatrix& operator=(const QRMatrix&) = delete;


public:

    // Constructors

        //- No default construct
        QRMatrix() = delete;

        //- Construct with a matrix and perform QR decomposition
        explicit QRMatrix
        (
            const modes mode,
            const outputs output,
            const bool pivoting,
            MatrixType& A
        );

        //- Construct with a const matrix and perform QR decomposition
        explicit QRMatrix
        (
            const MatrixType& A,
            const modes mode,
            const outputs output = outputs::BOTH_QR,
            const bool pivoting = false
        );


    // Member Functions

    // Access

        //- Return const reference to Q
        const MatrixType& Q() const noexcept
        {
            return Q_;
        }

        //- Return reference to Q
        MatrixType& Q() noexcept
        {
            return Q_;
        }

        //- Return const reference to R
        const MatrixType& R() const noexcept
        {
            return R_;
        }

        //- Return reference to R
        MatrixType& R() noexcept
        {
            return R_;
        }

        //- Return const reference to p
        const labelList& p() const noexcept
        {
            return p_;
        }

        //- Create and return the permutation matrix
        inline SMatrix P() const;


    // Linear system solution

        //- Solve the linear system with the given source
        //- and return the solution in the argument x
        void solve
        (
            List<cmptType>& x,
            const UList<cmptType>& source
        ) const;

        //- Solve the linear system with the given source
        //- and return the solution in the argument x
        template<class Addr>
        void solve
        (
            List<cmptType>& x,
            const IndirectListBase<cmptType, Addr>& source
        ) const;

        //- Solve the linear system with the given source
        //- and return the solution
        tmp<Field<cmptType>> solve
        (
            const UList<cmptType>& source
        ) const;

        //- Solve the linear system with the given source
        //- and return the solution
        template<class Addr>
        tmp<Field<cmptType>> solve
        (
            const IndirectListBase<cmptType, Addr>& source
        ) const;

        //- Solve a row-echelon-form linear system (Ax = b)
        //- starting from the bottom by back substitution
        //  A = R: Non-singular upper-triangular square matrix (m-by-m)
        //  b: Source (m-by-p)
        //  x: Solution (m-by-p)
        RMatrix solve
        (
            const RMatrix& b
        );

        //- Return the inverse of (Q*R), solving x = (Q*R).inv()*source
        SMatrix inv() const;
};

namespace MatrixTools
{

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Moore-Penrose inverse of singular/non-singular square/rectangular
//- scalar/complex matrices (KPP:p. 9834; KP:p. 648)
//  The tolerance to ensure the R1 matrix full-rank is set to 1e-5
//  by (TA; mentioned in (KPP:p. 9832)) in contrast to 1e-13 (KPP:p. 9834).
template<class MatrixType>
MatrixType pinv
(
    const MatrixType& A,
    scalar tol = 1e-5
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace MatrixTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "QRMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "QRMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
