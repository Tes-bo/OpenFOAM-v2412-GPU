/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "GeometricSphericalTensorField.H"
#include "sphericalTensorFieldField.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Cmpt, template<class> class PatchField, class GeoMesh>
void Foam::zip
(
    GeometricField<SphericalTensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Cmpt, PatchField, GeoMesh>& ii
)
{
    Foam::zip(result.primitiveFieldRef(), ii.primitiveField());

    Foam::zip(result.boundaryFieldRef(), ii.boundaryField());
}


template<class Cmpt, template<class> class PatchField, class GeoMesh>
void Foam::unzip
(
    const GeometricField<SphericalTensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Cmpt, PatchField, GeoMesh>& ii
)
{
    Foam::unzip(input.primitiveField(), ii.primitiveFieldRef());

    Foam::unzip(input.boundaryField(), ii.boundaryFieldRef());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, sphericalTensor, tr, transform)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph, transform)
UNARY_FUNCTION(scalar, sphericalTensor, det, pow3)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, inv, inv)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, '|', divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
