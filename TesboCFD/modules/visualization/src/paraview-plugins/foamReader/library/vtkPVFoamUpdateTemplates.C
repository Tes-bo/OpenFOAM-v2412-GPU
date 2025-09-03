/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoamUpdateTemplates_C
#define vtkPVFoamUpdateTemplates_C

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class patchType, class meshType>
void Foam::vtkPVFoam::updateInfoFields
(
    vtkDataArraySelection* select,
    const IOobjectList& objects
)
{
    DebugInfo
        << "<beg> updateInfoFields <"
        << meshType::Mesh::typeName
        << "> [volMeshPtr=" << (volMeshPtr_ ? "set" : "null") << "]"
        << nl;

    // Add geometric fields (volume/area) to GUI
    addToSelection<GeometricField<scalar, patchType, meshType>>
    (
        select,
        objects
    );
    addToSelection<GeometricField<vector, patchType, meshType>>
    (
        select,
        objects
    );
    addToSelection<GeometricField<sphericalTensor, patchType, meshType>>
    (
        select,
        objects
    );
    addToSelection<GeometricField<symmTensor, patchType, meshType>>
    (
        select,
        objects
    );
    addToSelection<GeometricField<tensor, patchType, meshType>>
    (
        select,
        objects
    );

    // Add dimensioned fields (volume/area) to GUI
    addToSelection<DimensionedField<scalar, meshType>>
    (
        select,
        objects
    );
    addToSelection<DimensionedField<vector, meshType>>
    (
        select,
        objects
    );
    addToSelection<DimensionedField<sphericalTensor, meshType>>
    (
        select,
        objects
    );
    addToSelection<DimensionedField<symmTensor, meshType>>
    (
        select,
        objects
    );
    addToSelection<DimensionedField<tensor, meshType>>
    (
        select,
        objects
    );

    DebugInfo << "<end> updateInfoFields" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
