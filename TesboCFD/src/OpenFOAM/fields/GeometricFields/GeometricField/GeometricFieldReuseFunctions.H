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

#ifndef Foam_GeometricFieldReuseFunctions_H
#define Foam_GeometricFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Detail
{

//- True if tmp GeometricField has a reusable pointer.
//
//  With GeometricField::debug on, check if patches have reusable types
template<class Type, template<class> class PatchField, class GeoMesh>
inline bool reusable
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    const bool ok = tfld.movable();

    if (ok && GeometricField<Type, PatchField, GeoMesh>::debug)
    {
        for (const auto& p : tfld().boundaryField())
        {
            if
            (
                !polyPatch::constraintType(p.patch().type())
             && !isA<typename PatchField<Type>::Calculated>(p)
            )
            {
                WarningInFunction
                    << "Attempt to reuse temporary with non-reusable BC "
                    << p.type() << endl;

                return false;
            }
        }
    }

    return ok;
}

} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// One-parameter versions

template
<
    class TypeR,
    class Type1,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpGeometricField
{
    //- Pass-through to New GeometricField
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const GeometricField<Type1, PatchField, GeoMesh>& f1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        return GeometricField<TypeR, PatchField, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }

    //- Dissimilar types: return new field
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& f1 = tf1();

        return GeometricField<TypeR, PatchField, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
struct reuseTmpGeometricField<TypeR, TypeR, PatchField, GeoMesh>
{
    //- Allow optional copy assignment of the initial content
    //- for identical input and output types
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tf1,
        const word& name,
        const dimensionSet& dimensions,
        const bool initCopy = false
    )
    {
        if (Detail::reusable(tf1))
        {
            auto& f1 = tf1.constCast();

            f1.rename(name);
            f1.dimensions().reset(dimensions);
            return tf1;
        }

        const auto& f1 = tf1();

        auto tresult = GeometricField<TypeR, PatchField, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );

        if (initCopy)
        {
            tresult.ref() == tf1();
        }

        return tresult;
    }
};


//- This global function forwards to reuseTmpGeometricField::New
template<class TypeR, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField<TypeR, PatchField, GeoMesh>
> New
(
    const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tf1,
    const word& name,
    const dimensionSet& dimensions,
    const bool initCopy = false
)
{
    return reuseTmpGeometricField<TypeR, TypeR, PatchField, GeoMesh>::New
    (
        tf1,
        name,
        dimensions,
        initCopy
    );
}


template
<
    class TypeR,
    class Type1,
    class Type12,
    class Type2,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpTmpGeometricField
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,
        const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& f1 = tf1();

        return GeometricField<TypeR, PatchField, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template
<
    class TypeR,
    class Type1,
    class Type12,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpTmpGeometricField
<
    TypeR, Type1, Type12, TypeR, PatchField, GeoMesh
>
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (Detail::reusable(tf2))
        {
            auto& f2 = tf2.constCast();

            f2.rename(name);
            f2.dimensions().reset(dimensions);
            return tf2;
        }

        const auto& f1 = tf1();

        return GeometricField<TypeR, PatchField, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template
<
    class TypeR,
    class Type2,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpTmpGeometricField
<
    TypeR, TypeR, TypeR, Type2, PatchField, GeoMesh
>
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tf1,
        const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (Detail::reusable(tf1))
        {
            auto& f1 = tf1.constCast();

            f1.rename(name);
            f1.dimensions().reset(dimensions);
            return tf1;
        }

        const auto& f1 = tf1();

        return GeometricField<TypeR, PatchField, GeoMesh>::New
        (
            name,
            f1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
struct reuseTmpTmpGeometricField
<
    TypeR, TypeR, TypeR, TypeR, PatchField, GeoMesh
>
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tf1,
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (Detail::reusable(tf1))
        {
            auto& f1 = tf1.constCast();

            f1.rename(name);
            f1.dimensions().reset(dimensions);
            return tf1;
        }
        if (Detail::reusable(tf2))
        {
            auto& f2 = tf2.constCast();

            f2.rename(name);
            f2.dimensions().reset(dimensions);
            return tf2;
        }

        const auto& f1 = tf1();

        return GeometricField<TypeR, PatchField, GeoMesh>::New
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
