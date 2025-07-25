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

#ifndef Foam_surfFieldsFwd_H
#define Foam_surfFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class GeoMesh> class DimensionedField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Face Fields

class surfGeoMesh;

// Typedefs
typedef
    DimensionedField<label, surfGeoMesh>
    surfLabelField;

typedef
    DimensionedField<scalar, surfGeoMesh>
    surfScalarField;

typedef
    DimensionedField<vector, surfGeoMesh>
    surfVectorField;

typedef
    DimensionedField<sphericalTensor, surfGeoMesh>
    surfSphericalTensorField;

typedef
    DimensionedField<symmTensor, surfGeoMesh>
    surfSymmTensorField;

typedef
    DimensionedField<tensor, surfGeoMesh>
    surfTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Point Fields

class surfPointGeoMesh;

// Typedefs
typedef
    DimensionedField<label, surfPointGeoMesh>
    surfPointLabelField;

typedef
    DimensionedField<scalar, surfPointGeoMesh>
    surfPointScalarField;

typedef
    DimensionedField<vector, surfPointGeoMesh>
    surfPointVectorField;

typedef
    DimensionedField<sphericalTensor, surfPointGeoMesh>
    surfPointSphericalTensorField;

typedef
    DimensionedField<symmTensor, surfPointGeoMesh>
    surfPointSymmTensorField;

typedef
    DimensionedField<tensor, surfPointGeoMesh>
    surfPointTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
