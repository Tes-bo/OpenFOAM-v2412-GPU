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

\*---------------------------------------------------------------------------*/

#ifndef Foam_FieldReuseFunctions_H
#define Foam_FieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// One-parameter versions

template<class TypeR, class Type1>
struct reuseTmp
{
    //- Pass-through to tmp New
    static tmp<Field<TypeR>> New(const Field<Type1>& f1)
    {
        return tmp<Field<TypeR>>::New(f1.size());
    }

    //- Dissimilar types: just use size
    static tmp<Field<TypeR>> New(const tmp<Field<Type1>>& tf1)
    {
        return tmp<Field<TypeR>>::New(tf1().size());
    }
};


template<class TypeR>
struct reuseTmp<TypeR, TypeR>
{
    //- Identical input and return types:
    //- allow optional copy assignment of the initial content
    static tmp<Field<TypeR>> New
    (
        const tmp<Field<TypeR>>& tf1,
        const bool initCopy = false
    )
    {
        if (tf1.movable())
        {
            return tf1;
        }

        auto rtf = tmp<Field<TypeR>>::New(tf1().size());

        if (initCopy)
        {
            rtf.ref() = tf1();
        }

        return rtf;
    }
};


//- This global function forwards to reuseTmp::New
template<class TypeR> tmp<Field<TypeR>> New
(
    const tmp<Field<TypeR>>& tf1,
    const bool initCopy = false
)
{
    return reuseTmp<TypeR, TypeR>::New(tf1, initCopy);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Two-parameter versions

template<class TypeR, class Type1, class Type12, class Type2>
struct reuseTmpTmp
{
    //- Dissimilar types: just use size
    static tmp<Field<TypeR>> New
    (
        const tmp<Field<Type1>>& tf1,
        const tmp<Field<Type2>>& tf2
    )
    {
        return tmp<Field<TypeR>>::New(tf1().size());
    }
};


template<class TypeR, class Type1, class Type12>
struct reuseTmpTmp<TypeR, Type1, Type12, TypeR>
{
    //- Second input has return type
    static tmp<Field<TypeR>> New
    (
        const tmp<Field<Type1>>& tf1,
        const tmp<Field<TypeR>>& tf2
    )
    {
        if (tf2.movable())
        {
            return tf2;
        }

        return tmp<Field<TypeR>>::New(tf1().size());
    }
};


template<class TypeR, class Type2>
struct reuseTmpTmp<TypeR, TypeR, TypeR, Type2>
{
    //- First input has return type
    static tmp<Field<TypeR>> New
    (
        const tmp<Field<TypeR>>& tf1,
        const tmp<Field<Type2>>& tf2
    )
    {
        if (tf1.movable())
        {
            return tf1;
        }

        return tmp<Field<TypeR>>::New(tf1().size());
    }
};


template<class TypeR>
struct reuseTmpTmp<TypeR, TypeR, TypeR, TypeR>
{
    //- Both inputs have return type
    static tmp<Field<TypeR>> New
    (
        const tmp<Field<TypeR>>& tf1,
        const tmp<Field<TypeR>>& tf2
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

        return tmp<Field<TypeR>>::New(tf1().size());
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
