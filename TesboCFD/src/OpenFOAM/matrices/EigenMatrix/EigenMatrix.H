/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::EigenMatrix

Description
    EigenMatrix (i.e. eigendecomposition or spectral decomposition) decomposes
    a diagonalisable nonsymmetric real square matrix into its canonical form,
    whereby the matrix is represented in terms of its eigenvalues and
    eigenvectors.

    The eigenvalue equation (i.e. eigenvalue problem) is:

    \f[
        A v = \lambda v
    \f]

    where
    \vartable
      A       | a diagonalisable square matrix of dimension m-by-m
      v       | a (non-zero) vector of dimension m (right eigenvector)
      \lambda | a scalar corresponding to v (eigenvalue)
    \endvartable


    If \c A is symmetric, the following relation is satisfied:

    \f[
        A = v*D*v^T
    \f]

    where
    \vartable
      D | diagonal real eigenvalue matrix
      v | orthogonal eigenvector matrix
    \endvartable

    If \c A is not symmetric, \c D becomes a block diagonal matrix wherein
    the real eigenvalues are present on the diagonal within 1-by-1 blocks, and
    complex eigenvalues within 2-by-2 blocks, i.e. \f$\lambda + i \mu\f$ with
    \f$[\lambda, \mu; -\mu, \lambda]\f$.

    The columns of \c v represent eigenvectors corresponding to eigenvalues,
    satisfying the eigenvalue equation. Even though eigenvalues of a matrix
    are unique, eigenvectors of the matrix are not. For the same eigenvalue,
    the corresponding eigenvector can be real or complex with non-unique
    entries. In addition, the validity of the equation \f$A = v*D*v^T\f$
    depends on the condition number of \c v, which can be ill-conditioned,
    or singular for invalidated equations.

    References:
    \verbatim
        OpenFOAM-compatible implementation:
            Passalacqua, A., Heylmun, J., Icardi, M.,
            Madadi, E., Bachant, P., & Hu, X. (2019).
            OpenQBMM 5.0.1 for OpenFOAM 7, Zenodo.
            DOI:10.5281/zenodo.3471804

        Implementations for the functions:
        'tridiagonaliseSymmMatrix', 'symmTridiagonalQL',
        'Hessenberg' and 'realSchur' (based on ALGOL-procedure:tred2):
            Wilkinson, J. H., & Reinsch, C. (1971).
            In Bauer, F. L. & Householder A. S. (Eds.),
            Handbook for Automatic Computation: Volume II: Linear Algebra.
            (Vol. 186), Springer-Verlag Berlin Heidelberg.
            DOI: 10.1007/978-3-642-86940-2

        Explanations on how real eigenvectors
        can be unpacked into complex domain:
            Moler, C. (1998).
            Re: Eigenvectors.
            Retrieved from https://bit.ly/3ao4Wat

        TNT/JAMA implementation:
            Pozo, R. (1997).
            Template Numerical Toolkit for linear algebra:
            High performance programming with C++
            and the Standard Template Library.
            The International Journal of Supercomputer Applications
            and High Performance Computing, 11(3), 251-263.
            DOI:10.1177/109434209701100307

            (No particular order) Hicklin, J., Moler, C., Webb, P.,
            Boisvert, R. F., Miller, B., Pozo, R., & Remington, K. (2012).
            JAMA: A Java Matrix Package 1.0.3.
            Retrived from https://math.nist.gov/javanumerics/jama/
    \endverbatim

Note
    - This implementation is an integration of the \c OpenQBMM \c eigenSolver
    class (2019) without any changes to its internal mechanisms. Therefore, no
    differences between EigenMatrix and \c eigenSolver (2019) classes should be
    expected in terms of input-process-output operations.
    - The \c OpenQBMM \c eigenSolver class derives almost completely from the
    \c TNT/JAMA implementation, a public-domain library developed by \c NIST
    and \c MathWorks from 1998 to 2012, available at
    http://math.nist.gov/tnt/index.html (Retrieved June 6, 2020). Their
    implementation was based upon \c EISPACK.
    - The \c tridiagonaliseSymmMatrix, \c symmTridiagonalQL, \c Hessenberg and
    \c realSchur methods are based on the \c Algol procedures \c tred2 by
    Bowdler, Martin, Reinsch, and Wilkinson, Handbook for Auto. Comp.,
    Vol. II-Linear Algebra, and the corresponding \c FORTRAN subroutine
    in \c EISPACK.

See also
    Test-EigenMatrix.C

SourceFiles
    EigenMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_EigenMatrix_H
#define Foam_EigenMatrix_H

#include "scalarMatrices.H"
#include "DiagonalMatrix.H"
#include "SquareMatrix.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class EigenMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class cmptType>
class EigenMatrix
{
    static_assert
    (
        std::is_floating_point<cmptType>::value,
        "EigenMatrix operates only with scalar base type."
    );

    // Private Data

        //- Number of rows and columns in input matrix
        const label n_;

        //- Real eigenvalues or real part of complex eigenvalues
        DiagonalMatrix<cmptType> EValsRe_;

        //- Zero-matrix for real eigenvalues
        //- or imaginary part of complex eigenvalues
        DiagonalMatrix<cmptType> EValsIm_;

        //- Right eigenvectors matrix where each column is
        //- a right eigenvector that corresponds to an eigenvalue
        SquareMatrix<cmptType> EVecs_;

        //- Copy of nonsymmetric input matrix evolving to eigenvectors matrix
        SquareMatrix<cmptType> H_;


    // Private Member Functions

        //- Householder transform of a symmetric matrix to tri-diagonal form
        void tridiagonaliseSymmMatrix();

        //- Symmetric tri-diagonal QL algorithm
        void symmTridiagonalQL();

        //- Reduce non-symmetric matrix to upper-Hessenberg form
        void Hessenberg();

        //- Reduce matrix from Hessenberg to real Schur form
        void realSchur();


public:

    // Generated Methods

        //- No default construct
        EigenMatrix() = delete;

        //- No copy construct
        EigenMatrix(const EigenMatrix&) = delete;

        //- No copy assignment
        EigenMatrix& operator=(const EigenMatrix&) = delete;


    // Constructors

        //- Construct from a SquareMatrix<cmptType>
        explicit EigenMatrix(const SquareMatrix<cmptType>& A);

        //- Construct from a SquareMatrix<cmptType> and symmetry flag
        //  Does not perform symmetric check
        EigenMatrix(const SquareMatrix<cmptType>& A, bool symmetric);


    // Access

        //- Return real eigenvalues or real part of complex eigenvalues
        const DiagonalMatrix<cmptType>& EValsRe() const noexcept
        {
            return EValsRe_;
        }

        //- Return zero-matrix for real eigenvalues
        //- or imaginary part of complex eigenvalues
        const DiagonalMatrix<cmptType>& EValsIm() const noexcept
        {
            return EValsIm_;
        }

        //- Return right eigenvectors matrix where each column is
        //- a right eigenvector that corresponds to an eigenvalue
        const SquareMatrix<cmptType>& EVecs() const noexcept
        {
            return EVecs_;
        }

        //- Return right eigenvectors in unpacked form
        const SquareMatrix<complex> complexEVecs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EigenMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
