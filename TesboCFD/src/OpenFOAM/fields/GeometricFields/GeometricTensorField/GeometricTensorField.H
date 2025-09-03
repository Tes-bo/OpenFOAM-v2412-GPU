/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

InClass
    Foam::GeometricTensorField

Description
    Tensor specific part of the implementation of GeometricField.

SourceFiles
    GeometricTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricTensorField_H
#define Foam_GeometricTensorField_H

#include "GeometricSphericalTensorField.H"
#include "GeometricSymmTensorField.H"
#include "DimensionedTensorField.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Zip together tensor field from components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void zip
(
    GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xx,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xy,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xz,
    const GeometricField<Cmpt, PatchField, GeoMesh>& yx,
    const GeometricField<Cmpt, PatchField, GeoMesh>& yy,
    const GeometricField<Cmpt, PatchField, GeoMesh>& yz,
    const GeometricField<Cmpt, PatchField, GeoMesh>& zx,
    const GeometricField<Cmpt, PatchField, GeoMesh>& zy,
    const GeometricField<Cmpt, PatchField, GeoMesh>& zz
);

//- Unzip tensor field into components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzip
(
    const GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Cmpt, PatchField, GeoMesh>& xx,
    GeometricField<Cmpt, PatchField, GeoMesh>& xy,
    GeometricField<Cmpt, PatchField, GeoMesh>& xz,
    GeometricField<Cmpt, PatchField, GeoMesh>& yx,
    GeometricField<Cmpt, PatchField, GeoMesh>& yy,
    GeometricField<Cmpt, PatchField, GeoMesh>& yz,
    GeometricField<Cmpt, PatchField, GeoMesh>& zx,
    GeometricField<Cmpt, PatchField, GeoMesh>& zy,
    GeometricField<Cmpt, PatchField, GeoMesh>& zz
);


//- Zip together tensor field from row components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void zipRows
(
    GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);

//- Zip together tensor field from column components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void zipCols
(
    GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);

//- Extract tensor field rows
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipRows
(
    const GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);

//- Extract tensor field columns
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipCols
(
    const GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);


//- Extract a tensor field row (x,y,z) == (0,1,2)
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipRow
(
    const GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& input,
    const direction idx,  //!< vector::components
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result
);

//- Extract a tensor field column (x,y,z) == (0,1,2)
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipCol
(
    const GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& input,
    const direction idx,  //!< vector::components
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result
);

//- Extract a tensor field diagonal
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipDiag
(
    const GeometricField<Tensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(tensor, tensor, T, transform)
UNARY_FUNCTION(scalar, tensor, tr, transform)
UNARY_FUNCTION(sphericalTensor, tensor, sph, transform)
UNARY_FUNCTION(symmTensor, tensor, symm, transform)
UNARY_FUNCTION(symmTensor, tensor, twoSymm, transform)
UNARY_FUNCTION(symmTensor, tensor, devSymm, transform)
UNARY_FUNCTION(symmTensor, tensor, devTwoSymm, transform)
UNARY_FUNCTION(tensor, tensor, skew, transform)
UNARY_FUNCTION(tensor, tensor, dev, transform)
UNARY_FUNCTION(tensor, tensor, dev2, transform)
UNARY_FUNCTION(scalar, tensor, det, transform)
UNARY_FUNCTION(tensor, tensor, cof, cof)
UNARY_FUNCTION(tensor, tensor, inv, inv)

UNARY_FUNCTION(vector, symmTensor, eigenValues, transform)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors, sign)


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual, transform)
UNARY_OPERATOR(tensor, vector, *, hdual, transform)

BINARY_OPERATOR(vector, vector, tensor, /, '|', divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#ifdef NoRepository
    #include "GeometricTensorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
