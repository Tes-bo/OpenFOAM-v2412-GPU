/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Macros for creating pointPatchField types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointPatchFieldMacros_H
#define Foam_pointPatchFieldMacros_H

#include "fieldTypes.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#define addToPointPatchFieldRunTimeSelection\
(PatchTypeField, typePatchTypeField)\
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        patch                                                                 \
    );                                                                        \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        patchMapper                                                           \
    );                                                                        \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        dictionary                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// For non-templated (or typedef) patch fields.
//- Define a concrete pointPatchField type and add to run-time tables
//- Example, (pointPatchScalarField, calculatedPointPatchScalarField)
#define makePointPatchTypeField(PatchTypeField, typePatchTypeField)           \
                                                                              \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addToPointPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField);

//- Define a pointPatchField type for a specific field type
//- Example, (scalar, calculated)
#define makePointPatchFieldTypeName(fieldType, bcType)                        \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CAT4(bcType, PointPatch, CAPITALIZE(fieldType), Field),               \
        0                                                                     \
    );

//- Add runTime selection for pointPatchField type of a specific field type
//- Example, (scalar, calculated)
#define addPointPatchFieldTypeRunTime(fieldType, bcType)                      \
    addToPointPatchFieldRunTimeSelection                                      \
    (                                                                         \
        CAT3(pointPatch, CAPITALIZE(fieldType), Field),                       \
        CAT4(bcType, PointPatch, CAPITALIZE(fieldType), Field)                \
    );

// For templated patch fields
#define makePointPatchFieldType(fieldType, bcType)                            \
    makePointPatchFieldTypeName(fieldType, bcType)                            \
    addPointPatchFieldTypeRunTime(fieldType, bcType)

//- Declare a typedef for a pointPatchField.
//- Example, (scalar, calculated)
#define makePointPatchFieldTypedef(fieldType, bcType)                         \
    typedef                                                                   \
        bcType##PointPatchField<fieldType>                                    \
        CAT4(bcType, PointPatch, CAPITALIZE(fieldType), Field);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- For all standard field types:
//- define a pointPatchField type and add to run-time tables
#define makePointPatchFields(bcType)                                          \
                                                                              \
    makePointPatchFieldType(scalar, bcType);                                  \
    makePointPatchFieldType(vector, bcType);                                  \
    makePointPatchFieldType(sphericalTensor, bcType);                         \
    makePointPatchFieldType(symmTensor, bcType);                              \
    makePointPatchFieldType(tensor, bcType);


//- For all standard field types: define a pointPatchField type
//- Example, (calculated)
#define makePointPatchFieldsTypeName(bcType)                                  \
                                                                              \
    makePointPatchFieldTypeName(scalar, bcType);                              \
    makePointPatchFieldTypeName(vector, bcType);                              \
    makePointPatchFieldTypeName(sphericalTensor, bcType);                     \
    makePointPatchFieldTypeName(symmTensor, bcType);                          \
    makePointPatchFieldTypeName(tensor, bcType);


//- For all standard field types: define typedefs for the pointPatchField
//- Example, (calculated)
#define makePointPatchFieldTypedefs(bcType)                                   \
                                                                              \
    makePointPatchFieldTypedef(scalar, bcType);                               \
    makePointPatchFieldTypedef(vector, bcType);                               \
    makePointPatchFieldTypedef(sphericalTensor, bcType);                      \
    makePointPatchFieldTypedef(symmTensor, bcType);                           \
    makePointPatchFieldTypedef(tensor, bcType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
