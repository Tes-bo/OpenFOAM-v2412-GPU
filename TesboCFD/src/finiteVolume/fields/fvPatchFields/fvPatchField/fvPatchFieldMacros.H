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
    Macros for creating fvPatchField types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvPatchFieldMacros_H
#define Foam_fvPatchFieldMacros_H

#include "fieldTypes.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#undef  addToPatchFieldRunTimeSelection
#define addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)   \
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


// Use with caution
#undef  addRemovableToPatchFieldRunTimeSelection
#define addRemovableToPatchFieldRunTimeSelection\
(PatchTypeField, typePatchTypeField)                                          \
                                                                              \
    addRemovableToRunTimeSelectionTable                                       \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        patch                                                                 \
    );                                                                        \
    addRemovableToRunTimeSelectionTable                                       \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        patchMapper                                                           \
    );                                                                        \
    addRemovableToRunTimeSelectionTable                                       \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        dictionary                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// For non-templated (or typedef) patch fields.
//- Define a concrete fvPatchField type and add to run-time tables
//- Example, (fvPatchScalarField, calculatedFvPatchScalarField)
#define makePatchTypeField(PatchTypeField, typePatchTypeField)                \
                                                                              \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField);

// For non-templated patch fields - use with caution
#define makeRemovablePatchTypeField(PatchTypeField, typePatchTypeField)       \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addRemovableToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


//- Define a fvPatchField type for a specific field type
//- Example, (scalar, calculated)
#define makePatchFieldTypeName(fieldType, bcType)                             \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CAT4(bcType, FvPatch, CAPITALIZE(fieldType), Field),                  \
        0                                                                     \
    );

//- Add runTime selection for fvPatchField type of a specific field type
//- Example, (scalar, calculated)
#define addPatchFieldTypeRunTime(fieldType, bcType)                           \
    addToPatchFieldRunTimeSelection                                           \
    (                                                                         \
        CAT3(fvPatch, CAPITALIZE(fieldType), Field),                          \
        CAT4(bcType, FvPatch, CAPITALIZE(fieldType), Field)                   \
    );

// For templated patch fields
#define makePatchFieldType(fieldType, bcType)                                 \
    makePatchFieldTypeName(fieldType, bcType)                                 \
    addPatchFieldTypeRunTime(fieldType, bcType)

//- Declare a typedef for a fvPatchField. Example, (scalar, calculated)
#define makePatchTypeFieldTypedef(fieldType, bcType)                          \
    typedef                                                                   \
        bcType##FvPatchField<fieldType>                                       \
        CAT4(bcType, FvPatch, CAPITALIZE(fieldType), Field);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- For all standard field types:
//- define a faPatchField type and add to run-time tables
#define makePatchFields(bcType)                                               \
                                                                              \
    makePatchFieldType(scalar, bcType);                                       \
    makePatchFieldType(vector, bcType);                                       \
    makePatchFieldType(sphericalTensor, bcType);                              \
    makePatchFieldType(symmTensor, bcType);                                   \
    makePatchFieldType(tensor, bcType);


//- For all standard field types: define a faPatchField type.
//- Example, (calculated)
#define makePatchFieldTypeNames(bcType)                                       \
                                                                              \
    makePatchFieldTypeName(scalar, bcType);                                   \
    makePatchFieldTypeName(vector, bcType);                                   \
    makePatchFieldTypeName(sphericalTensor, bcType);                          \
    makePatchFieldTypeName(symmTensor, bcType);                               \
    makePatchFieldTypeName(tensor, bcType);


//- For all standard field types: define typedefs for the faPatchField
//- Example, (calculated)
#define makePatchTypeFieldTypedefs(bcType)                                    \
                                                                              \
    makePatchTypeFieldTypedef(scalar, bcType);                                \
    makePatchTypeFieldTypedef(vector, bcType);                                \
    makePatchTypeFieldTypedef(sphericalTensor, bcType);                       \
    makePatchTypeFieldTypedef(symmTensor, bcType);                            \
    makePatchTypeFieldTypedef(tensor, bcType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compatibility: OpenFOAM-v2212 and earlier

#define makeTemplatePatchTypeField(fldType, bcType) \
    makePatchFieldType(fldType, bcType)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
