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
    Forwards and collection of common sliced edge field types

\*---------------------------------------------------------------------------*/

#ifndef Foam_slicedEdgeFieldsFwd_H
#define Foam_slicedEdgeFieldsFwd_H

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
class edgeMesh;
template<class Type> class faePatchField;
template<class Type> class slicedFaePatchField;


//- A (volume) surface field slice for a given type
template<class Type>
using EdgeSliceField =
SlicedGeometricField<Type, faePatchField, slicedFaePatchField, edgeMesh>;


// Typedefs

typedef
    SlicedGeometricField
    <scalar, faePatchField, slicedFaePatchField, edgeMesh>
    slicedEdgeScalarField;

typedef
    SlicedGeometricField
    <vector, faePatchField, slicedFaePatchField, edgeMesh>
    slicedEdgeVectorField;

typedef
    SlicedGeometricField
    <sphericalTensor, faePatchField, slicedFaePatchField, edgeMesh>
    slicedEdgeSphericalTensorField;

typedef
    SlicedGeometricField
    <symmTensor, faePatchField, slicedFaePatchField, edgeMesh>
    slicedEdgeSymmTensorField;

typedef
    SlicedGeometricField
    <tensor, faePatchField, slicedFaePatchField, edgeMesh>
    slicedEdgeTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
