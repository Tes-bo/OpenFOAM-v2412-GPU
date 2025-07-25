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

InClass
    Foam::surfaceFields

Description

SourceFiles
    surfaceFields.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceFieldsFwd_H
#define Foam_surfaceFieldsFwd_H

#include "fieldTypes.H"
#include "GeometricFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class surfaceMesh;
template<class Type> class fvsPatchField;


//- A (volume) surface field for a given type
template<class Type>
using SurfaceField = GeometricField<Type, fvsPatchField, surfaceMesh>;


// Typedefs

typedef
    GeometricField<scalar, fvsPatchField, surfaceMesh>
    surfaceScalarField;

typedef
    GeometricField<vector, fvsPatchField, surfaceMesh>
    surfaceVectorField;

typedef
    GeometricField<sphericalTensor, fvsPatchField, surfaceMesh>
    surfaceSphericalTensorField;

typedef
    GeometricField<symmTensor, fvsPatchField, surfaceMesh>
    surfaceSymmTensorField;

typedef
    GeometricField<tensor, fvsPatchField, surfaceMesh>
    surfaceTensorField;


/*---------------------------------------------------------------------------*\
                       Namespace fieldTypes Declaration
\*---------------------------------------------------------------------------*/

namespace fieldTypes
{
    //- Standard surface field types (scalar, vector, tensor, etc)
    extern const wordList surface;

} // End namespace fieldTypes


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
