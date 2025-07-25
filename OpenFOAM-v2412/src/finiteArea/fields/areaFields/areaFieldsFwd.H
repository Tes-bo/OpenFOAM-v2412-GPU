/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

Description
    Forwards and collection of common area field types

SourceFiles
    areaFields.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_areaFieldsFwd_H
#define Foam_areaFieldsFwd_H

#include "fieldTypes.H"
#include "GeometricFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class areaMesh;
template<class Type> class faPatchField;


//- An area field for a given type
template<class Type>
using AreaField = GeometricField<Type, faPatchField, areaMesh>;

//- An area boundary field for a given Type
template<class Type>
using AreaBoundaryField = GeometricBoundaryField<Type, faPatchField, areaMesh>;

//- An area internal field for a given type
template<class Type>
using AreaInternalField = DimensionedField<Type, areaMesh>;


// Typedefs

typedef
    GeometricField<scalar, faPatchField, areaMesh>
    areaScalarField;

typedef
    GeometricField<vector, faPatchField, areaMesh>
    areaVectorField;

typedef
    GeometricField<sphericalTensor, faPatchField, areaMesh>
    areaSphericalTensorField;

typedef
    GeometricField<symmTensor, faPatchField, areaMesh>
    areaSymmTensorField;

typedef
    GeometricField<tensor, faPatchField, areaMesh>
    areaTensorField;


/*---------------------------------------------------------------------------*\
                       Namespace fieldTypes Declaration
\*---------------------------------------------------------------------------*/

namespace fieldTypes
{
    //- Standard area field types (scalar, vector, tensor, etc)
    extern const wordList area;

    //- Standard dimensioned field types (scalar, vector, tensor, etc)
    extern const wordList area_internal;

} // End namespace fieldTypes


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
