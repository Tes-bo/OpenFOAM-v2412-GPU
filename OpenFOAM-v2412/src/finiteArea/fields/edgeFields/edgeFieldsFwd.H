/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

InNamespace
    Foam

Description
    Forwards for edge field types

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeFieldsFwd_H
#define Foam_edgeFieldsFwd_H

#include "fieldTypes.H"
#include "GeometricFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class edgeMesh;
template<class Type> class faePatchField;


/// Future? (TBD)
/// //- An (area) edge field for a given type
/// template<class Type>
/// using EdgeField = GeometricField<Type, faePatchField, edgeMesh>;


// Typedefs

typedef
    GeometricField<scalar, faePatchField, edgeMesh>
    edgeScalarField;

typedef
    GeometricField<vector, faePatchField, edgeMesh>
    edgeVectorField;

typedef
    GeometricField<sphericalTensor, faePatchField, edgeMesh>
    edgeSphericalTensorField;

typedef
    GeometricField<symmTensor, faePatchField, edgeMesh>
    edgeSymmTensorField;

typedef
    GeometricField<tensor, faePatchField, edgeMesh>
    edgeTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
