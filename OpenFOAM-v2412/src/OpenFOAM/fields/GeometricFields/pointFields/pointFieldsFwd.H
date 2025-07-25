/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::pointFields

Description
    Forwards and collection of common point field types

SourceFiles
    pointFields.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointFieldsFwd_H
#define Foam_pointFieldsFwd_H

#include "fieldTypes.H"
#include "GeometricFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class pointMesh;
template<class Type> class pointPatchField;


//- A point field for a given type
template<class Type>
using PointField = GeometricField<Type, pointPatchField, pointMesh>;

//- A point boundary field for a given Type
template<class Type>
using PointBoundaryField =
    GeometricBoundaryField<Type, pointPatchField, pointMesh>;

//- A point internal field for a given type
template<class Type>
using PointInternalField = DimensionedField<Type, pointMesh>;


// Typedefs

typedef
    GeometricField<scalar, pointPatchField, pointMesh>
    pointScalarField;

typedef
    GeometricField<vector, pointPatchField, pointMesh>
    pointVectorField;

typedef
    GeometricField<sphericalTensor, pointPatchField, pointMesh>
    pointSphericalTensorField;

typedef
    GeometricField<symmTensor, pointPatchField, pointMesh>
    pointSymmTensorField;

typedef
    GeometricField<tensor, pointPatchField, pointMesh>
    pointTensorField;


/*---------------------------------------------------------------------------*\
                       Namespace fieldTypes Declaration
\*---------------------------------------------------------------------------*/

namespace fieldTypes
{
    //- Standard point field types (scalar, vector, tensor, etc)
    extern const wordList point;

} // End namespace fieldTypes


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
