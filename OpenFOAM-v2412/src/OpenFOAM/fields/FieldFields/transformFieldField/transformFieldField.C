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

\*---------------------------------------------------------------------------*/

#include "transformFieldField.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
void Foam::transform
(
    FieldField<Field, Type>& result,
    const tensor& rot,
    const FieldField<Field, Type>& fld
)
{
    forAll(result, i)
    {
        transform(result[i], rot, fld[i]);
    }
}


template<template<class> class Field, class Type>
void Foam::transform
(
    FieldField<Field, Type>& result,
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
)
{
    forAll(result, i)
    {
        transform(result[i], rot[i], fld[i]);
    }
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::transform
(
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    transform(tranf(), rot, fld);
    return tranf;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::transform
(
    const FieldField<Field, tensor>& rot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    transform(tresult(), rot, tresult());
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::transform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tresult
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    transform(tresult(), trot(), fld);
    trot.clear();
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::transform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    transform(tresult(), trot(), tresult());
    trot.clear();
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::transform
(
    const tensor& rot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tresult
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    transform(tresult(), rot, fld);
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::transform
(
    const tensor& rot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    transform(tresult(), rot, tresult());
    return tresult;
}


template<template<class> class Field, class Type>
void Foam::invTransform
(
    FieldField<Field, Type>& result,
    const tensor& rot,
    const FieldField<Field, Type>& fld
)
{
    forAll(result, i)
    {
        invTransform(result[i], rot, fld[i]);
    }
}


template<template<class> class Field, class Type>
void Foam::invTransform
(
    FieldField<Field, Type>& result,
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
)
{
    forAll(result, i)
    {
        invTransform(result[i], rot[i], fld[i]);
    }
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::invTransform
(
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    invTransform(tranf(), rot, fld);
    return tranf;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::invTransform
(
    const FieldField<Field, tensor>& rot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    invTransform(tresult(), rot, tresult());
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tresult
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    invTransform(tresult(), trot(), fld);
    trot.clear();
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    invTransform(tresult(), trot(), tresult());
    trot.clear();
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::invTransform
(
    const tensor& rot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tresult
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    invTransform(tresult(), rot, fld);
    return tresult;
}


template<template<class> class Field, class Type>
Foam::tmp<Foam::FieldField<Field, Type>>
Foam::invTransform
(
    const tensor& rot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    invTransform(tresult(), rot, tresult());
    return tresult;
}


// ************************************************************************* //
