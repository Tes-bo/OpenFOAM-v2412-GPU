/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "fileFieldSelection.H"
#include "objectRegistry.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "fvPatchField.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "GeometricField.H"
#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fileFieldSelection::addFromFile
(
    const IOobjectList& objects,
    DynamicList<fieldInfo>& set
) const
{
    for (const fieldInfo& fi : *this)
    {
        const wordList names(objects.sortedNames<Type>(fi.name()));

        if (!names.empty())
        {
            fi.found(true);
        }
        for (const word& name : names)
        {
            set.emplace_back(wordRe(name));
        }
    }
}


template<template<class> class PatchType, class MeshType>
void Foam::functionObjects::fileFieldSelection::addGeoFieldTypes
(
    const IOobjectList& objects,
    DynamicList<fieldInfo>& set
) const
{
    #undef  doLocalCode
    #define doLocalCode(DataType)                                             \
    addFromFile<GeometricField<DataType, PatchType, MeshType>>(objects, set);

    doLocalCode(scalar);
    doLocalCode(vector);
    doLocalCode(sphericalTensor);
    doLocalCode(symmTensor);
    doLocalCode(tensor);
    #undef doLocalCode
}


void Foam::functionObjects::fileFieldSelection::addInternalFieldTypes
(
    const IOobjectList& objects,
    DynamicList<fieldInfo>& set
) const
{
    #undef  doLocalCode
    #define doLocalCode(DataType)                                             \
    addFromFile<DimensionedField<DataType, volMesh>>(objects, set);

    doLocalCode(scalar);
    doLocalCode(vector);
    doLocalCode(sphericalTensor);
    doLocalCode(symmTensor);
    doLocalCode(tensor);
    #undef doLocalCode
}


void Foam::functionObjects::fileFieldSelection::addUniformFieldTypes
(
    const IOobjectList& objects,
    DynamicList<fieldInfo>& set
) const
{
    #undef  doLocalCode
    #define doLocalCode(DataType)                                             \
    addFromFile<UniformDimensionedField<DataType>>(objects, set);

    doLocalCode(scalar);
    doLocalCode(vector);
    doLocalCode(sphericalTensor);
    doLocalCode(symmTensor);
    doLocalCode(tensor);
    #undef doLocalCode
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fileFieldSelection::fileFieldSelection
(
    const objectRegistry& obr,
    const bool includeComponents
)
:
    fieldSelection(obr, includeComponents)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fileFieldSelection::updateSelection()
{
    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);
    const IOobjectList objects(mesh, mesh.time().timeName());

    List<fieldInfo> oldSet(std::move(selection_));

    DynamicList<fieldInfo> newSelection(oldSet.size());

    // Geometric fields
    addGeoFieldTypes<fvPatchField, volMesh>(objects, newSelection);
    addGeoFieldTypes<fvsPatchField, surfaceMesh>(objects, newSelection);
    addGeoFieldTypes<pointPatchField, pointMesh>(objects, newSelection);

    // Internal fields
    addInternalFieldTypes(objects, newSelection);

    // Uniform fields
    addUniformFieldTypes(objects, newSelection);

    selection_.transfer(newSelection);

    (void)fieldSelection::checkSelection();

    return selection_ != oldSet;
}


// ************************************************************************* //
