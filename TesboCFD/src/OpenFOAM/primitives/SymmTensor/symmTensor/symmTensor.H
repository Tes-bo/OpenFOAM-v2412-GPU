/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

Typedef
    Foam::symmTensor

Description
    SymmTensor of scalars, i.e. SymmTensor<scalar>.

    Analytical functions for the computation of real eigenvalues and
    real eigenvectors from a given symmTensor.

See also
    Test-SymmTensor.C

SourceFiles
    symmTensor.C

\*---------------------------------------------------------------------------*/

#ifndef symmTensor_H
#define symmTensor_H

#include "SymmTensor.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef SymmTensor<scalar> symmTensor;

typedef Tensor<scalar> tensor;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return real ascending-order eigenvalues of a given symmTensor
//  \param T symmTensor
//
//  \return vector real eigenvalues
vector eigenValues(const symmTensor& T);


//- Return a real eigenvector corresponding to
//- a given real eigenvalue of a given symmTensor
//  \param T symmTensor
//  \param eVal real eigenvalue
//  \param standardBasis1 symmTensor orthogonal component 1
//  \param standardBasis2 symmTensor orthogonal component 2
//
//  \return vector real eigenvector
vector eigenVector
(
    const symmTensor& T,
    const scalar eVal,
    const vector& standardBasis1,
    const vector& standardBasis2
);


//- Return real eigenvectors corresponding to
//- given real eigenvalues of a given symmTensor
//  \param T symmTensor
//  \param eVals real eigenvalues
//
//  \return tensor real eigenvectors, each row is an eigenvector
tensor eigenVectors
(
    const symmTensor& T,
    const vector& eVals
);


//- Return real eigenvectors of a given symmTensor by computing
//- the real eigenvalues of the tensor in the background
//  \param T symmTensor
//
//  \return tensor real eigenvectors, each row is an eigenvector
tensor eigenVectors(const symmTensor& T);


//- Return inverse of a given symmTensor, and fall back
//- to pseudo-inverse if the symmTensor is singular
//  \param st symmTensor
//
//  \return symmTensor inverse of symmTensor
symmTensor pinv(const symmTensor& st);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
