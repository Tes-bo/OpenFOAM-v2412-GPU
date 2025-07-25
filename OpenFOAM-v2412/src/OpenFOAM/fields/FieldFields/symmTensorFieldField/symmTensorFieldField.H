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
    Foam::symmTensorFieldField

Description
    Specialisation of FieldField\<T\> for symmTensor.

SourceFiles
    symmTensorFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_symmTensorFieldField_H
#define Foam_symmTensorFieldField_H

#include "FieldField.H"
#include "symmTensor.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Zip together symmTensor field field from components
template<template<class> class Field, class Cmpt>
void zip
(
    FieldField<Field, SymmTensor<Cmpt>>& result,
    const FieldField<Field, Cmpt>& xx,
    const FieldField<Field, Cmpt>& xy,
    const FieldField<Field, Cmpt>& xz,
    const FieldField<Field, Cmpt>& yy,
    const FieldField<Field, Cmpt>& yz,
    const FieldField<Field, Cmpt>& zz
);

//- Unzip symmTensor field field into components
template<template<class> class Field, class Cmpt>
void unzip
(
    const FieldField<Field, SymmTensor<Cmpt>>& input,
    FieldField<Field, Cmpt>& xx,
    FieldField<Field, Cmpt>& xy,
    FieldField<Field, Cmpt>& xz,
    FieldField<Field, Cmpt>& yy,
    FieldField<Field, Cmpt>& yz,
    FieldField<Field, Cmpt>& zz
);

//- Zip together symmTensor field field from row components
template<template<class> class Field, class Cmpt>
void zipRows
(
    FieldField<Field, SymmTensor<Cmpt>>& result,
    const FieldField<Field, Vector<Cmpt>>& x,
    const FieldField<Field, Vector<Cmpt>>& y,
    const FieldField<Field, Vector<Cmpt>>& z
);

//- Zip together symmTensor field from column components
template<template<class> class Field, class Cmpt>
void zipCols
(
    FieldField<Field, SymmTensor<Cmpt>>& result,
    const FieldField<Field, Vector<Cmpt>>& x,
    const FieldField<Field, Vector<Cmpt>>& y,
    const FieldField<Field, Vector<Cmpt>>& z
);


//- Extract symmTensor field field rows
template<template<class> class Field, class Cmpt>
void unzipRows
(
    const FieldField<Field, SymmTensor<Cmpt>>& input,
    FieldField<Field, Vector<Cmpt>>& x,
    FieldField<Field, Vector<Cmpt>>& y,
    FieldField<Field, Vector<Cmpt>>& z
);

//- Extract symmTensor field field columns
template<template<class> class Field, class Cmpt>
void unzipCols
(
    const FieldField<Field, SymmTensor<Cmpt>>& input,
    FieldField<Field, Vector<Cmpt>>& x,
    FieldField<Field, Vector<Cmpt>>& y,
    FieldField<Field, Vector<Cmpt>>& z
);


//- Extract a symmTensor field field row (x,y,z) == (0,1,2)
template<template<class> class Field, class Cmpt>
void unzipRow
(
    const FieldField<Field, SymmTensor<Cmpt>>& input,
    const direction idx,  //!< vector::components
    FieldField<Field, Vector<Cmpt>>& result
);

//- Extract a symmTensor field field column (x,y,z) == (0,1,2)
template<template<class> class Field, class Cmpt>
void unzipCol
(
    const FieldField<Field, SymmTensor<Cmpt>>& input,
    const direction idx,  //!< vector::components
    FieldField<Field, Vector<Cmpt>>& result
);

//- Extract a symmTensor field field diagonal
template<template<class> class Field, class Cmpt>
void unzipDiag
(
    const FieldField<Field, SymmTensor<Cmpt>>& input,
    FieldField<Field, Vector<Cmpt>>& result
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)
UNARY_FUNCTION(symmTensor, symmTensor, innerSqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, devSymm)
UNARY_FUNCTION(symmTensor, symmTensor, devTwoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)
UNARY_FUNCTION(symmTensor, symmTensor, inv)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#ifdef NoRepository
    #include "symmTensorFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
