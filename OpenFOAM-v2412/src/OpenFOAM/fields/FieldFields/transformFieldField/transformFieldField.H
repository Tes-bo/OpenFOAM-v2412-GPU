/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    transformFieldField
    Spatial transformation functions for FieldField.

SourceFiles
    transformFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_transformFieldField_H
#define Foam_transformFieldField_H

#include "transform.H"
#include "tensorFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// transform()

template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>& result,
    const tensor& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>& result,
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>& rot,
    const tmp<FieldField<Field, Type>>& tfld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const tmp<FieldField<Field, Type>>& tfld
);


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor& rot,
    const tmp<FieldField<Field, Type>>& tfld
);


// invTransform()

template<template<class> class Field, class Type>
void invTransform
(
    FieldField<Field, Type>& result,
    const tensor& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
void invTransform
(
    FieldField<Field, Type>& result,
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const FieldField<Field, tensor>& rot,
    const tmp<FieldField<Field, Type>>& tfld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const tmp<FieldField<Field, Type>>& tfld
);


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tensor& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tensor& rot,
    const tmp<FieldField<Field, Type>>& tfld
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "transformFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
