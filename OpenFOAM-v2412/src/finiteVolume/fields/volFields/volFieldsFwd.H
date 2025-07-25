/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Forwards and collection of common volume field types

SourceFiles
    volFields.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_volFieldsFwd_H
#define Foam_volFieldsFwd_H

#include "fieldTypes.H"
#include "GeometricFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class volMesh;
template<class Type> class fvPatchField;


//- A volume field for a given type
template<class Type>
using VolumeField = GeometricField<Type, fvPatchField, volMesh>;

//- A volume boundary field for a given Type
template<class Type>
using VolumeBoundaryField = GeometricBoundaryField<Type, fvPatchField, volMesh>;

//- A volume internal field for a given Type
template<class Type>
using VolumeInternalField = DimensionedField<Type, volMesh>;


// Typedefs

typedef
    GeometricField<scalar, fvPatchField, volMesh>
    volScalarField;

typedef
    GeometricField<vector, fvPatchField, volMesh>
    volVectorField;

typedef
    GeometricField<sphericalTensor, fvPatchField, volMesh>
    volSphericalTensorField;

typedef
    GeometricField<symmTensor, fvPatchField, volMesh>
    volSymmTensorField;

typedef
    GeometricField<tensor, fvPatchField, volMesh>
    volTensorField;


/*---------------------------------------------------------------------------*\
                       Namespace fieldTypes Declaration
\*---------------------------------------------------------------------------*/

namespace fieldTypes
{
    //- Standard dimensioned field types (scalar, vector, tensor, etc)
    extern const wordList internal;

    //- Standard volume field types (scalar, vector, tensor, etc)
    extern const wordList volume;

} // End namespace fieldTypes


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
