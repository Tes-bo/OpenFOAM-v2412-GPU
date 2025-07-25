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

InClass
    Foam::transformGeometricField

Description
    Spatial transformation functions for GeometricField

SourceFiles
    transformGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_transformGeometricField_H
#define Foam_transformGeometricField_H

#include "transform.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// transform()

template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const dimensionedTensor& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
);



// invTransform()

template<class Type, template<class> class PatchField, class GeoMesh>
void invTransform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
void invTransform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
);


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const dimensionedTensor& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "transformGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
