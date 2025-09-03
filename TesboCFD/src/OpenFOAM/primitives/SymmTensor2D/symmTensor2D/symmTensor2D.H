/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

Typedef
    Foam::symmTensor2D

Description
    SymmTensor2D of scalars, i.e. SymmTensor2D<scalar>.

    Analytical functions for the computation of real eigenvalues
    and real eigenvectors from a given symmTensor2D.

    Reference:
    \verbatim
        2-by-2 eigenvalue algorithm (tag:K):
            Kronenburg, M. J. (2015).
            A method for fast diagonalization of
            a 2x2 or 3x3 real symmetric matrix.
            arXiv:1306.6291.
    \endverbatim

See also
    Test-SymmTensor2D.C

SourceFiles
    symmTensor2D.C

\*---------------------------------------------------------------------------*/

#ifndef symmTensor2D_H
#define symmTensor2D_H

#include "SymmTensor2D.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef SymmTensor2D<scalar> symmTensor2D;

typedef Tensor2D<scalar> tensor2D;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return real arbitrary-order eigenvalues of a given symmTensor2D
//  \param T symmTensor2D
//
//  \return vector2D real eigenvalues
vector2D eigenValues(const symmTensor2D& T);


//- Return a real eigenvector corresponding to
//- a given real eigenvalue of a given symmTensor2D
//  \param T symmTensor2D
//  \param eVal real eigenvalue
//  \param standardBasis symmTensor2D orthogonal component, e.g. vector2D(1, 0)
//
//  \return vector2D real eigenvector
vector2D eigenVector
(
    const symmTensor2D& T,
    const scalar eVal,
    const vector2D& standardBasis
);


//- Return real eigenvectors corresponding to
//- given real eigenvalues of a given symmTensor2D
//  \param T symmTensor2D
//  \param eVals real eigenvalues
//
//  \return tensor2D real eigenvectors, each row is an eigenvector
tensor2D eigenVectors(const symmTensor2D& T, const vector2D& eVals);


//- Return real eigenvectors of a given symmTensor2D by computing
//- the real eigenvalues of the symmTensor2D in the background
//  \param T symmTensor2D
//
//  \return tensor2D real eigenvectors, each row is an eigenvector
tensor2D eigenVectors(const symmTensor2D& T);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
