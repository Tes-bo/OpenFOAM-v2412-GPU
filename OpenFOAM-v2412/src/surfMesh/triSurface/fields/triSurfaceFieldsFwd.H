/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#ifndef Foam_triSurfaceFieldsFwd_H
#define Foam_triSurfaceFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class GeoMesh> class DimensionedField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Face Fields

class triSurfaceGeoMesh;

// Typedefs
typedef
    DimensionedField<label, triSurfaceGeoMesh>
    triSurfaceLabelField;

typedef
    DimensionedField<scalar, triSurfaceGeoMesh>
    triSurfaceScalarField;

typedef
    DimensionedField<vector, triSurfaceGeoMesh>
    triSurfaceVectorField;

typedef
    DimensionedField<sphericalTensor, triSurfaceGeoMesh>
    triSurfaceSphericalTensorField;

typedef
    DimensionedField<symmTensor, triSurfaceGeoMesh>
    triSurfaceSymmTensorField;

typedef
    DimensionedField<tensor, triSurfaceGeoMesh>
    triSurfaceTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Point Fields

class triSurfacePointGeoMesh;

typedef
    DimensionedField<label, triSurfacePointGeoMesh>
    triSurfacePointLabelField;

typedef
    DimensionedField<scalar, triSurfacePointGeoMesh>
    triSurfacePointScalarField;

typedef
    DimensionedField<vector, triSurfacePointGeoMesh>
    triSurfacePointVectorField;

typedef
    DimensionedField<sphericalTensor, triSurfacePointGeoMesh>
    triSurfacePointSphericalTensorField;

typedef
    DimensionedField<symmTensor, triSurfacePointGeoMesh>
    triSurfacePointSymmTensorField;

typedef
    DimensionedField<tensor, triSurfacePointGeoMesh>
    triSurfacePointTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
