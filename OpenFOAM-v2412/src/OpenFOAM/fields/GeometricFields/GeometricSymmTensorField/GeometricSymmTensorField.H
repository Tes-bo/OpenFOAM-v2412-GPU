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
    Foam::GeometricSymmTensorField

Description
    SymmTensor specific part of the implementation of GeometricField.

SourceFiles
    GeometricSymmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricSymmTensorField_H
#define Foam_GeometricSymmTensorField_H

#include "GeometricField.H"
#include "DimensionedSymmTensorField.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Zip together symmTensor field from components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void zip
(
    GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xx,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xy,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xz,
    const GeometricField<Cmpt, PatchField, GeoMesh>& yy,
    const GeometricField<Cmpt, PatchField, GeoMesh>& yz,
    const GeometricField<Cmpt, PatchField, GeoMesh>& zz
);

//- Unzip symmTensor field into components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzip
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Cmpt, PatchField, GeoMesh>& xx,
    GeometricField<Cmpt, PatchField, GeoMesh>& xy,
    GeometricField<Cmpt, PatchField, GeoMesh>& xz,
    GeometricField<Cmpt, PatchField, GeoMesh>& yy,
    GeometricField<Cmpt, PatchField, GeoMesh>& yz,
    GeometricField<Cmpt, PatchField, GeoMesh>& zz
);


//- Zip together symmTensor field from row components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void zipRows
(
    GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);

//- Zip together symmTensor field from column components
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void zipCols
(
    GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);

//- Extract symmTensor field rows
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipRows
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);

//- Extract symmTensor field columns
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipCols
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& x,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& y,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& z
);


//- Extract a symmTensor field row (x,y,z) == (0,1,2)
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipRow
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    const direction idx,  //!< vector::components
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result
);

//- Extract a symmTensor field column (x,y,z) == (0,1,2)
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipCol
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    const direction idx,  //!< vector::components
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result
);


//- Extract a symmTensor field diagonal
template<class Cmpt, template<class> class PatchField, class GeoMesh>
void unzipDiag
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr, sqr)
UNARY_FUNCTION(symmTensor, symmTensor, innerSqr, sqr)

UNARY_FUNCTION(scalar, symmTensor, tr, transform)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph, transform)
UNARY_FUNCTION(symmTensor, symmTensor, symm, transform)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm, transform)
UNARY_FUNCTION(symmTensor, symmTensor, devSymm, transform)
UNARY_FUNCTION(symmTensor, symmTensor, devTwoSymm, transform)
UNARY_FUNCTION(symmTensor, symmTensor, dev, transform)
UNARY_FUNCTION(symmTensor, symmTensor, dev2, transform)
UNARY_FUNCTION(scalar, symmTensor, det, transform)
UNARY_FUNCTION(symmTensor, symmTensor, cof, cof)
UNARY_FUNCTION(symmTensor, symmTensor, inv, inv)


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual, transform)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, '&', dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, '&', dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#ifdef NoRepository
    #include "GeometricSymmTensorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
