/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Wikki Ltd
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
    Macros for creating faPatchField types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_faPatchFieldMacros_H
#define Foam_faPatchFieldMacros_H

#include "fieldTypes.H"
#include "faPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#define addToFaPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField) \
                                                                              \
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
//- Define a concrete faPatchField type and add to run-time tables
//- Example, (faPatchScalarField, calculatedFaPatchScalarField)
#define makeFaPatchTypeField(PatchTypeField, typePatchTypeField)              \
                                                                              \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addToFaPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField);

//- Define a faPatchField type for a specific field type
//- Example, (scalar, calculated)
#define makeFaPatchFieldTypeName(fieldType, bcType)                           \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CAT4(bcType, FaPatch, CAPITALIZE(fieldType), Field),                  \
        0                                                                     \
    );

//- Add runTime selection for faPatchField type of a specific field type
//- Example, (scalar, calculated)
#define addFaPatchFieldTypeRunTime(fieldType, bcType)                         \
    addToFaPatchFieldRunTimeSelection                                         \
    (                                                                         \
        CAT3(faPatch, CAPITALIZE(fieldType), Field),                          \
        CAT4(bcType, FaPatch, CAPITALIZE(fieldType), Field)                   \
    );

// For templated patch fields
#define makeFaPatchFieldType(fieldType, bcType)                               \
    makeFaPatchFieldTypeName(fieldType, bcType)                               \
    addFaPatchFieldTypeRunTime(fieldType, bcType)                             \

//- Declare a typedef for a faPatchField. Example, (scalar, calculated)
#define makeFaPatchFieldTypedef(fieldType, bcType)                            \
    typedef                                                                   \
        bcType##FaPatchField<fieldType>                                       \
        CAT4(bcType, FaPatch, CAPITALIZE(fieldType), Field);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- For all standard field types:
//- define a faPatchField type and add to run-time tables
#define makeFaPatchFields(bcType)                                             \
                                                                              \
    makeFaPatchFieldType(scalar, bcType);                                     \
    makeFaPatchFieldType(vector, bcType);                                     \
    makeFaPatchFieldType(sphericalTensor, bcType);                            \
    makeFaPatchFieldType(symmTensor, bcType);                                 \
    makeFaPatchFieldType(tensor, bcType);


//- For all standard field types: define a faPatchField type.
//- Example, (calculated)
#define makeFaPatchFieldsTypeName(bcType)                                     \
                                                                              \
    makeFaPatchFieldTypeName(scalar, bcType);                                 \
    makeFaPatchFieldTypeName(vector, bcType);                                 \
    makeFaPatchFieldTypeName(sphericalTensor, bcType);                        \
    makeFaPatchFieldTypeName(symmTensor, bcType);                             \
    makeFaPatchFieldTypeName(tensor, bcType);


//- For all standard field types: define typedefs for the faPatchField.
//- Example, (calculated)
#define makeFaPatchTypeFieldTypedefs(bcType)                                  \
                                                                              \
    makeFaPatchFieldTypedef(scalar, bcType);                                  \
    makeFaPatchFieldTypedef(vector, bcType);                                  \
    makeFaPatchFieldTypedef(sphericalTensor, bcType);                         \
    makeFaPatchFieldTypedef(symmTensor, bcType);                              \
    makeFaPatchFieldTypedef(tensor, bcType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
