/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

Global
    Foam::ReadFields

Description
    Field reading functions for post-processing utilities

SourceFiles
    ReadFields.C
    ReadFieldsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ReadFields_H
#define Foam_ReadFields_H

#include "PtrList.H"
#include "wordList.H"
#include "DynamicList.H"
#include "GeometricField.H"
#include "HashSet.H"
#include "LIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regIOobject;
class IOobjectList;
class objectRegistry;

//- Read Geometric fields of templated type.
//  \return sorted names of fields read.
//  \note All processors guaranteed to read fields in same order.
template<class Type, template<class> class PatchField, class GeoMesh>
wordList ReadFields
(
    const typename GeoMesh::Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
    const bool syncPar = true,
    const bool readOldTime = false
);


//- Read fields of the templated type.
//  \return sorted names of fields read.
//  \note All processors guaranteed to read fields in same order.
template<class GeoField, class Mesh>
wordList ReadFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar = true
);

//- Read non-mesh fields (uniformDimensionedField like 'g').
//  \return sorted names of fields read.
//  \note All processors guaranteed to read fields in same order.
template<class GeoField>
wordList ReadFields
(
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar = true
);

//- Read all GeometricFields of the templated type.
//  \param fieldsCache is an objectRegistry of all stored fields
template<class GeoField>
static void ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    objectRegistry& fieldsCache
);

//- Read all GeometricFields of the templated type.
//  \param fieldsCache is the objectRegistry name where fields are stored
template<class GeoField>
static void ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    const word& registryName = "fieldsCache"
);


//- Read the selected GeometricFields of the templated type
//- and store on the objectRegistry.
//  Returns a list of field pointers for later cleanup
template<class GeoFieldType, class NameMatchPredicate>
void readFields
(
    const typename GeoFieldType::Mesh& mesh,
    const IOobjectList& objects,
    //! Restrict to fields with matching names
    const NameMatchPredicate& selectedFields,
    //! [out] List of field pointers for later cleanup
    DynamicList<regIOobject*>& storedObjects
);

//- Read the selected UniformDimensionedFields of the templated type
//- and store on the objectRegistry.
//  Returns a list of field pointers for later cleanup
template<class UniformFieldType, class NameMatchPredicate>
void readUniformFields
(
    const IOobjectList& objects,
    //! Restrict to fields with matching names
    const NameMatchPredicate& selectedFields,
    //! [out] List of field pointers for later cleanup
    DynamicList<regIOobject*>& storedObjects,
    const bool syncPar = true
);


// Housekeeping

//- Read the selected GeometricFields of the templated type
//- and store on the objectRegistry.
//  \deprecated(2023-07) - prefer the DynamicList version
//  Returns a stack of field pointers for later cleanup
template<class GeoFieldType, class NameMatchPredicate>
FOAM_DEPRECATED_FOR(2023-07, "DynamicList version")
void readFields
(
    const typename GeoFieldType::Mesh& mesh,
    const IOobjectList& objects,
    const NameMatchPredicate& selectedFields,
    LIFOStack<regIOobject*>& storedObjects
);


//- Read the selected UniformDimensionedFields of the templated type
//- and store on the objectRegistry.
//  \deprecated(2023-07) - prefer the DynamicList version
//  Returns a stack of field pointers for later cleanup
template<class UniformFieldType, class NameMatchPredicate>
FOAM_DEPRECATED_FOR(2023-07, "DynamicList version")
void readUniformFields
(
    const IOobjectList& objects,
    const NameMatchPredicate& selectedFields,
    LIFOStack<regIOobject*>& storedObjects,
    const bool syncPar = true
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReadFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
