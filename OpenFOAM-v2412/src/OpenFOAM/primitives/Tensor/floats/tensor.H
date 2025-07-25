/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::tensor

Description
    Tensor of scalars, i.e. Tensor\<scalar\>.

    Analytical functions for the computation of complex eigenvalues and
    complex eigenvectors from a given tensor.

See also
    Test-Tensor.C

SourceFiles
    floatTensors.C
    tensor.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tensor_H
#define Foam_tensor_H

#include "Tensor.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//! \class Foam::floatTensor
//! \brief A Tensor of values with float precision
typedef Tensor<float> floatTensor;

//! \class Foam::doubleTensor
//! \brief A Tensor of values with double precision
typedef Tensor<double> doubleTensor;

// With float or double precision (depending on compilation)
typedef Tensor<scalar> tensor;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return complex eigenvalues of a given tensor
//  \param T tensor
//
//  \return Vector<complex> eigenvalues
Vector<complex> eigenValues(const tensor& T);


//- Return a complex eigenvector corresponding to
//- a given complex eigenvalue of a given tensor
//  \param T tensor
//  \param eVal complex eigenvalue
//  \param standardBasis1 tensor orthogonal component 1
//  \param standardBasis2 tensor orthogonal component 2
//
//  \return Vector<complex> eigenvector
Vector<complex> eigenVector
(
    const tensor& T,
    const complex& eVal,
    const Vector<complex>& standardBasis1,
    const Vector<complex>& standardBasis2
);


//- Return complex eigenvectors corresponding to
//- given complex eigenvalues of a given tensor
//  \param T tensor
//  \param eVals complex eigenvalues
//
//  \return Tensor<complex> eigenvectors, each row is an eigenvector
Tensor<complex> eigenVectors
(
    const tensor& T,
    const Vector<complex>& eVals
);


//- Return complex eigenvectors of a given tensor by computing
//- the complex eigenvalues of the tensor in the background
//  \param T tensor
//
//  \return Tensor<complex> complex eigenvectors, each row is an eigenvector
Tensor<complex> eigenVectors(const tensor& T);


//- Return inverse of a given tensor, and fall back
//- into pseudo-inverse if the tensor is singular
//  \param t tensor
//
//  \return tensor inverse of tensor
tensor pinv(const tensor& t);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
