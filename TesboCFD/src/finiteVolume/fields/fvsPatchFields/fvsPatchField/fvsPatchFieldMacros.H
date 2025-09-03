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
    Macros for creating fvsPatchField types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvsPatchFieldMacros_H
#define Foam_fvsPatchFieldMacros_H

#include "fieldTypes.H"
#include "fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#define addToFvsPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)\
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
//- Define a concrete fvsPatchField type and add to run-time tables
//- Example, (fvsPatchScalarField, calculatedFvsPatchScalarField)
#define makeFvsPatchTypeField(PatchTypeField, typePatchTypeField)             \
                                                                              \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addToFvsPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField);


//- Define a fvsPatchField type for a specific field type
//- Example, (scalar, calculated)
#define makeFvsPatchFieldTypeName(fieldType, bcType)                          \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CAT4(bcType, FvsPatch, CAPITALIZE(fieldType), Field),                 \
        0                                                                     \
    );

//- Add runTime selection for fvsPatchField type of a specific field type
//- Example, (scalar, calculated)
#define addFvsPatchFieldTypeRunTime(fieldType, bcType)                        \
    addToFvsPatchFieldRunTimeSelection                                        \
    (                                                                         \
        CAT3(fvsPatch, CAPITALIZE(fieldType), Field),                         \
        CAT4(bcType, FvsPatch, CAPITALIZE(fieldType), Field)                  \
    );

// For templated patch fields
#define makeFvsPatchFieldType(fieldType, bcType)                              \
    makeFvsPatchFieldTypeName(fieldType, bcType)                              \
    addFvsPatchFieldTypeRunTime(fieldType, bcType)

//- Declare a typedef for a fvsPatchField. Example, (scalar, calculated)
#define makeFvsPatchFieldTypedef(fieldType, bcType)                           \
    typedef                                                                   \
        bcType##FvsPatchField<fieldType>                                      \
        CAT4(bcType, FvsPatch, CAPITALIZE(fieldType), Field);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- For all standard field types:
//- define a fvsPatchField type and add to run-time tables
#define makeFvsPatchFields(bcType)                                            \
                                                                              \
    makeFvsPatchFieldType(scalar, bcType);                                    \
    makeFvsPatchFieldType(vector, bcType);                                    \
    makeFvsPatchFieldType(sphericalTensor, bcType);                           \
    makeFvsPatchFieldType(symmTensor, bcType);                                \
    makeFvsPatchFieldType(tensor, bcType);


//- For all standard field types: define a fvsPatchField type
//- Example, (calculated)
#define makeFvsPatchFieldsTypeName(bcType)                                    \
                                                                              \
    makeFvsPatchFieldTypeName(scalar, bcType);                                \
    makeFvsPatchFieldTypeName(vector, bcType);                                \
    makeFvsPatchFieldTypeName(sphericalTensor, bcType);                       \
    makeFvsPatchFieldTypeName(symmTensor, bcType);                            \
    makeFvsPatchFieldTypeName(tensor, bcType);


//- For all standard field types: define typedefs for the fvsPatchField
//- Example, (calculated)
#define makeFvsPatchTypeFieldTypedefs(bcType)                                 \
                                                                              \
    makeFvsPatchFieldTypedef(scalar, bcType);                                 \
    makeFvsPatchFieldTypedef(vector, bcType);                                 \
    makeFvsPatchFieldTypedef(sphericalTensor, bcType);                        \
    makeFvsPatchFieldTypedef(symmTensor, bcType);                             \
    makeFvsPatchFieldTypedef(tensor, bcType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
