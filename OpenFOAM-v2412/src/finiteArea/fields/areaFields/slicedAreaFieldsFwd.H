/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Forwards and collection of common sliced area field types

\*---------------------------------------------------------------------------*/

#ifndef Foam_slicedAreaFieldsFwd_H
#define Foam_slicedAreaFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
class SlicedGeometricField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class areaMesh;
template<class Type> class faPatchField;
template<class Type> class slicedFaPatchField;

//- An aree field slice for a given type
template<class Type>
using AreaSliceField =
SlicedGeometricField<Type, faPatchField, slicedFaPatchField, areaMesh>;


// Typedefs

typedef
    SlicedGeometricField
    <scalar, faPatchField, slicedFaPatchField, areaMesh>
    slicedAreaScalarField;

typedef
    SlicedGeometricField
    <vector, faPatchField, slicedFaPatchField, areaMesh>
    slicedAreaVectorField;

typedef
    SlicedGeometricField
    <sphericalTensor, faPatchField, slicedFaPatchField, areaMesh>
    slicedAreaSphericalTensorField;

typedef
    SlicedGeometricField
    <symmTensor, faPatchField, slicedFaPatchField, areaMesh>
    slicedAreaSymmTensorField;

typedef
    SlicedGeometricField
    <tensor, faPatchField, slicedFaPatchField, areaMesh>
    slicedAreaTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
