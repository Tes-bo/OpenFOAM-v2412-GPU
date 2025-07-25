/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "tensorFieldField.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.C"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<template<class> class Field, class Cmpt>
void Foam::zip
(
    FieldField<Field, Tensor<Cmpt>>& result,
    const FieldField<Field, Cmpt>& xx,
    const FieldField<Field, Cmpt>& xy,
    const FieldField<Field, Cmpt>& xz,
    const FieldField<Field, Cmpt>& yx,
    const FieldField<Field, Cmpt>& yy,
    const FieldField<Field, Cmpt>& yz,
    const FieldField<Field, Cmpt>& zx,
    const FieldField<Field, Cmpt>& zy,
    const FieldField<Field, Cmpt>& zz
)
{
    forAll(result, i)
    {
        Foam::zip
        (
            result[i],
            xx[i], xy[i], xz[i],
            yx[i], yy[i], yz[i],
            zx[i], zy[i], zz[i]
        );
    }
}


template<template<class> class Field, class Cmpt>
void Foam::unzip
(
    const FieldField<Field, Tensor<Cmpt>>& input,
    FieldField<Field, Cmpt>& xx,
    FieldField<Field, Cmpt>& xy,
    FieldField<Field, Cmpt>& xz,
    FieldField<Field, Cmpt>& yx,
    FieldField<Field, Cmpt>& yy,
    FieldField<Field, Cmpt>& yz,
    FieldField<Field, Cmpt>& zx,
    FieldField<Field, Cmpt>& zy,
    FieldField<Field, Cmpt>& zz
)
{
    forAll(input, i)
    {
        Foam::unzip
        (
            input[i],
            xx[i], xy[i], xz[i],
            yx[i], yy[i], yz[i],
            zx[i], zy[i], zz[i]
        );
    }
}


template<template<class> class Field, class Cmpt>
void Foam::zipRows
(
    FieldField<Field, Tensor<Cmpt>>& result,
    const FieldField<Field, Vector<Cmpt>>& x,
    const FieldField<Field, Vector<Cmpt>>& y,
    const FieldField<Field, Vector<Cmpt>>& z
)
{
    forAll(result, i)
    {
        Foam::zipRows(result[i], x[i], y[i], z[i]);
    }
}


template<template<class> class Field, class Cmpt>
void Foam::zipCols
(
    FieldField<Field, Tensor<Cmpt>>& result,
    const FieldField<Field, Vector<Cmpt>>& x,
    const FieldField<Field, Vector<Cmpt>>& y,
    const FieldField<Field, Vector<Cmpt>>& z
)
{
    forAll(result, i)
    {
        Foam::zipCols(result[i], x[i], y[i], z[i]);
    }
}


template<template<class> class Field, class Cmpt>
void Foam::unzipRows
(
    const FieldField<Field, Tensor<Cmpt>>& input,
    FieldField<Field, Vector<Cmpt>>& x,
    FieldField<Field, Vector<Cmpt>>& y,
    FieldField<Field, Vector<Cmpt>>& z
)
{
    forAll(input, i)
    {
        Foam::unzipRows(input[i], x[i], y[i], z[i]);
    }
}


template<template<class> class Field, class Cmpt>
void Foam::unzipCols
(
    const FieldField<Field, Tensor<Cmpt>>& input,
    FieldField<Field, Vector<Cmpt>>& x,
    FieldField<Field, Vector<Cmpt>>& y,
    FieldField<Field, Vector<Cmpt>>& z
)
{
    forAll(input, i)
    {
        Foam::unzipCols(input[i], x[i], y[i], z[i]);
    }
}


template<template<class> class Field, class Cmpt>
void Foam::unzipRow
(
    const FieldField<Field, Tensor<Cmpt>>& input,
    const direction idx,
    FieldField<Field, Vector<Cmpt>>& result
)
{
    forAll(input, i)
    {
        Foam::unzipRow(input[i], idx, result[i]);
    }
}


template<template<class> class Field, class Cmpt>
void Foam::unzipCol
(
    const FieldField<Field, Tensor<Cmpt>>& input,
    const direction idx,
    FieldField<Field, Vector<Cmpt>>& result
)
{
    forAll(input, i)
    {
        Foam::unzipCol(input[i], idx, result[i]);
    }
}


template<template<class> class Field, class Cmpt>
void Foam::unzipDiag
(
    const FieldField<Field, Tensor<Cmpt>>& input,
    FieldField<Field, Vector<Cmpt>>& result
)
{
    forAll(input, i)
    {
        Foam::unzipDiag(input[i], result[i]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(symmTensor, tensor, devSymm)
UNARY_FUNCTION(symmTensor, tensor, devTwoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)
UNARY_FUNCTION(tensor, tensor, inv)

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
