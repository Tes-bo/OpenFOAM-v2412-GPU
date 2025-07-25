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

Description
    Spatial transformation functions for FieldFields.

\*---------------------------------------------------------------------------*/

#include "transformGeometricField.H"
#include "transformField.H"
#include "transformFieldField.H"
#include "GeometricField.H"

// * * * * * * * * * * * * Transform Global Functions  * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::transform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    transform
    (
        result.primitiveFieldRef(),
        rot.value(),
        fld.primitiveField()
    );
    transform
    (
        result.boundaryFieldRef(),
        rot.value(),
        fld.boundaryField()
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::transform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    transform
    (
        result.primitiveFieldRef(),
        rot.primitiveField(),
        fld.primitiveField()
    );
    transform
    (
        result.boundaryFieldRef(),
        rot.boundaryField(),
        fld.boundaryField()
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = GeometricField<Type, PatchField, GeoMesh>::New
    (
        "transform(" + rot.name() + ',' + fld.name() + ')',
        fld.mesh(),
        fld.dimensions()
    );

    transform(tresult.ref(), rot, fld);

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = transform(rot, tfld());
    tfld.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = transform(trot(), fld);
    trot.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = transform(trot(), tfld());
    tfld.clear();
    trot.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::transform
(
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = GeometricField<Type, PatchField, GeoMesh>::New
    (
        "transform(" + rot.name() + ',' + fld.name() + ')',
        fld.mesh(),
        fld.dimensions()
    );

    transform(tresult.ref(), rot, fld);

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::transform
(
    const dimensionedTensor& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = transform(rot, tfld());
    tfld.clear();
    return tresult;
}


// * * * * * * * * * * * invTransform Global Functions * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::invTransform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    invTransform
    (
        result.primitiveFieldRef(),
        rot.value(),
        fld.primitiveField()
    );
    invTransform
    (
        result.boundaryFieldRef(),
        rot.value(),
        fld.boundaryField()
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::invTransform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    invTransform
    (
        result.primitiveFieldRef(),
        rot.primitiveField(),
        fld.primitiveField()
    );
    invTransform
    (
        result.boundaryFieldRef(),
        rot.boundaryField(),
        fld.boundaryField()
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = GeometricField<Type, PatchField, GeoMesh>::New
    (
        "invTransform(" + rot.name() + ',' + fld.name() + ')',
        fld.mesh(),
        fld.dimensions()
    );

    invTransform(tresult.ref(), rot, fld);

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = invTransform(rot, tfld());
    tfld.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = invTransform(trot(), fld);
    trot.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = invTransform(trot(), tfld());
    tfld.clear();
    trot.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::invTransform
(
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = GeometricField<Type, PatchField, GeoMesh>::New
    (
        "invTransform(" + rot.name() + ',' + fld.name() + ')',
        fld.mesh(),
        fld.dimensions()
    );

    invTransform(tresult.ref(), rot, fld);

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::invTransform
(
    const dimensionedTensor& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = invTransform(rot, tfld());
    tfld.clear();
    return tresult;
}


// ************************************************************************* //
