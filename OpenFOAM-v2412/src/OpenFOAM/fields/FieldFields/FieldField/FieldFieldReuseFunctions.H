/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#ifndef Foam_FieldFieldReuseFunctions_H
#define Foam_FieldFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// One-parameter versions

template<template<class> class Field, class TypeR, class Type1>
struct reuseTmpFieldField
{
    //- Pass-through to NewCalculatedType
    static tmp<FieldField<Field, TypeR>> New
    (
        const FieldField<Field, Type1>& f1
    )
    {
        return FieldField<Field, TypeR>::NewCalculatedType(f1);
    }

    //- Dissimilar types: just use size
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, Type1>>& tf1
    )
    {
        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR>
struct reuseTmpFieldField<Field, TypeR, TypeR>
{
    //- Identical input and return types:
    //- allow optional copy assignment of the initial content
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, TypeR>>& tf1,
        const bool initCopy = false
    )
    {
        if (tf1.movable())
        {
            return tf1;
        }

        auto tresult = FieldField<Field, TypeR>::NewCalculatedType(tf1());

        if (initCopy)
        {
            tresult.ref() = tf1();
        }

        return tresult;
    }
};


//- This global function forwards to reuseTmpFieldField::New
template<template<class> class Field, class TypeR>
tmp<FieldField<Field, TypeR>> New
(
    const tmp<FieldField<Field, TypeR>>& tf1,
    const bool initCopy = false
)
{
    return reuseTmpFieldField<Field, TypeR, TypeR>::New(tf1, initCopy);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Two-parameter versions

template
<
    template<class> class Field,
    class TypeR,
    class Type1,
    class Type12,
    class Type2
>
struct reuseTmpTmpFieldField
{
    //- Dissimilar types: just use size
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, Type1>>& tf1,
        const tmp<FieldField<Field, Type2>>& tf2
    )
    {
        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR, class Type1, class Type12>
struct reuseTmpTmpFieldField<Field, TypeR, Type1, Type12, TypeR>
{
    //- Second input has return type
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, Type1>>& tf1,
        const tmp<FieldField<Field, TypeR>>& tf2
    )
    {
        if (tf2.movable())
        {
            return tf2;
        }

        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR, class Type2>
struct reuseTmpTmpFieldField<Field, TypeR, TypeR, TypeR, Type2>
{
    //- First input has return type
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, TypeR>>& tf1,
        const tmp<FieldField<Field, Type2>>& tf2
    )
    {
        if (tf1.movable())
        {
            return tf1;
        }

        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR>
struct reuseTmpTmpFieldField<Field, TypeR, TypeR, TypeR, TypeR>
{
    //- Both inputs have return type
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, TypeR>>& tf1,
        const tmp<FieldField<Field, TypeR>>& tf2
    )
    {
        if (tf1.movable())
        {
            return tf1;
        }
        if (tf2.movable())
        {
            return tf2;
        }

        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
