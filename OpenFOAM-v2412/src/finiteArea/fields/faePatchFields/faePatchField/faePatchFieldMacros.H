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
    Macros for creating faePatchField types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_faePatchFieldMacros_H
#define Foam_faePatchFieldMacros_H

#include "fieldTypes.H"
#include "faePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#define addToFaePatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)\
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
//- Define a concrete faePatchField type and add to run-time tables
//- Example, (faePatchScalarField, calculatedFaePatchScalarField)
#define makeFaePatchTypeField(PatchTypeField, typePatchTypeField)             \
                                                                              \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addToFaePatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField);

//- Define a faePatchField type for a specific field type
//- Example, (scalar, calculated)
#define makeFaePatchFieldTypeName(fieldType, bcType)                          \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CAT4(bcType, FaePatch, CAPITALIZE(fieldType), Field),                 \
        0                                                                     \
    );

//- Add runTime selection for faePatchField type of a specific field type
//- Example, (scalar, calculated)
#define addFaePatchFieldTypeRunTime(fieldType, bcType)                        \
    addToFaePatchFieldRunTimeSelection                                        \
    (                                                                         \
        CAT3(faePatch, CAPITALIZE(fieldType), Field),                         \
        CAT4(bcType, FaePatch, CAPITALIZE(fieldType), Field)                  \
    );

// For templated patch fields
#define makeFaePatchFieldType(fieldType, bcType)                              \
    makeFaePatchFieldTypeName(fieldType, bcType)                              \
    addFaePatchFieldTypeRunTime(fieldType, bcType)

//- Declare a typedef for a faePatchField. Example, (scalar, calculated)
#define makeFaePatchFieldTypedef(fieldType, bcType)                           \
    typedef                                                                   \
        bcType##FaePatchField<fieldType>                                      \
        CAT4(bcType, FaePatch, CAPITALIZE(fieldType), Field);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- For all standard field types:
//- define a faePatchField type and add to run-time tables
#define makeFaePatchFields(bcType)                                            \
                                                                              \
    makeFaePatchFieldType(scalar, bcType);                                    \
    makeFaePatchFieldType(vector, bcType);                                    \
    makeFaePatchFieldType(sphericalTensor, bcType);                           \
    makeFaePatchFieldType(symmTensor, bcType);                                \
    makeFaePatchFieldType(tensor, bcType);


//- For all standard field types: define a faePatchField type
//- Example, (calculated)
#define makeFaePatchFieldsTypeName(bcType)                                    \
                                                                              \
    makeFaePatchFieldTypeName(scalar, bcType);                                \
    makeFaePatchFieldTypeName(vector, bcType);                                \
    makeFaePatchFieldTypeName(sphericalTensor, bcType);                       \
    makeFaePatchFieldTypeName(symmTensor, bcType);                            \
    makeFaePatchFieldTypeName(tensor, bcType);


//- For all standard field types: define typedefs for the faePatchField
//- Example, (calculated)
#define makeFaePatchTypeFieldTypedefs(bcType)                                 \
                                                                              \
    makeFaePatchFieldTypedef(scalar, bcType);                                 \
    makeFaePatchFieldTypedef(vector, bcType);                                 \
    makeFaePatchFieldTypedef(sphericalTensor, bcType);                        \
    makeFaePatchFieldTypedef(symmTensor, bcType);                             \
    makeFaePatchFieldTypedef(tensor, bcType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
