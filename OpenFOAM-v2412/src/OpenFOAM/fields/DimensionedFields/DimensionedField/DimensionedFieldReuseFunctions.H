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

#ifndef Foam_DimensionedFieldReuseFunctions_H
#define Foam_DimensionedFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// One-parameter versions

template<class TypeR, class Type1, class GeoMesh>
struct reuseTmpDimensionedField
{
    //- Pass-through to New DimensionedField
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const DimensionedField<Type1, GeoMesh>& f1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        return DimensionedField<TypeR, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }

    //- Dissimilar types: return new field
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh>>& tf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& f1 = tf1();

        return DimensionedField<TypeR, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class GeoMesh>
struct reuseTmpDimensionedField<TypeR, TypeR, GeoMesh>
{
    //- Allow optional copy assignment of the initial content
    //- for identical input and output types
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh>>& tf1,
        const word& name,
        const dimensionSet& dimensions,
        const bool initCopy = false
    )
    {
        if (tf1.movable())
        {
            auto& f1 = tf1.constCast();

            f1.rename(name);
            f1.dimensions().reset(dimensions);
            return tf1;
        }

        const auto& f1 = tf1();

        auto tresult = DimensionedField<TypeR, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );

        if (initCopy)
        {
            tresult.ref().field() = f1.field();
        }

        return tresult;
    }
};


//- Global function forwards to reuseTmpDimensionedField::New
template<class TypeR, class GeoMesh>
tmp<DimensionedField<TypeR, GeoMesh>> New
(
    const tmp<DimensionedField<TypeR, GeoMesh>>& tf1,
    const word& name,
    const dimensionSet& dimensions,
    const bool initCopy = false
)
{
    return reuseTmpDimensionedField<TypeR, TypeR, GeoMesh>::New
    (
        tf1,
        name,
        dimensions,
        initCopy
    );
}


template<class TypeR, class Type1, class Type12, class Type2, class GeoMesh>
struct reuseTmpTmpDimensionedField
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh>>& tf1,
        const tmp<DimensionedField<Type2, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& f1 = tf1();

        return DimensionedField<TypeR, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class Type1, class Type12, class GeoMesh>
struct reuseTmpTmpDimensionedField<TypeR, Type1, Type12, TypeR, GeoMesh>
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh>>& tf1,
        const tmp<DimensionedField<TypeR, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (tf2.movable())
        {
            auto& f2 = tf2.constCast();

            f2.rename(name);
            f2.dimensions().reset(dimensions);
            return tf2;
        }

        const auto& f1 = tf1();

        return DimensionedField<TypeR, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class Type2, class GeoMesh>
struct reuseTmpTmpDimensionedField<TypeR, TypeR, TypeR, Type2, GeoMesh>
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh>>& tf1,
        const tmp<DimensionedField<Type2, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (tf1.movable())
        {
            auto& f1 = tf1.constCast();

            f1.rename(name);
            f1.dimensions().reset(dimensions);
            return tf1;
        }

        const auto& f1 = tf1();

        return DimensionedField<TypeR, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class GeoMesh>
struct reuseTmpTmpDimensionedField<TypeR, TypeR, TypeR, TypeR, GeoMesh>
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh>>& tf1,
        const tmp<DimensionedField<TypeR, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (tf1.movable())
        {
            auto& f1 = tf1.constCast();

            f1.rename(name);
            f1.dimensions().reset(dimensions);
            return tf1;
        }
        if (tf2.movable())
        {
            auto& f2 = tf2.constCast();

            f2.rename(name);
            f2.dimensions().reset(dimensions);
            return tf2;
        }

        const auto& f1 = tf1();

        return DimensionedField<TypeR, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
