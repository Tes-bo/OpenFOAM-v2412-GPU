/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "polyMeshFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// Naming to shadow volScalarField::Internal etc.
// keep synchronized with finiteVolume volFields.C

template<>
const word DimensionedField<scalar, polyGeoMesh>::typeName
(
    "volScalarField::Internal"
);

template<>
const word DimensionedField<vector, polyGeoMesh>::typeName
(
    "volVectorField::Internal"
);

template<>
const word DimensionedField<sphericalTensor, polyGeoMesh>::typeName
(
    "volSphericalTensorField::Internal"
);

template<>
const word DimensionedField<symmTensor, polyGeoMesh>::typeName
(
    "volSymmTensorField::Internal"
);

template<>
const word DimensionedField<tensor, polyGeoMesh>::typeName
(
    "volTensorField::Internal"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
