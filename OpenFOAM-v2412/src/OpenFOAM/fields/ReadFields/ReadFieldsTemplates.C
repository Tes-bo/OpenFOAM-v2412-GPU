/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "ReadFields.H"
#include "HashSet.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::wordList Foam::ReadFields
(
    const typename GeoMesh::Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
    const bool syncPar,
    const bool readOldTime
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    // Names of GeoField objects, sorted order.
    const wordList fieldNames(objects.sortedNames<GeoField>(syncPar));

    // Construct the fields - reading in consistent (master) order.
    fields.resize_null(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Reading " << GeoField::typeName << ':';
        }
        Info<< ' ' << fieldName;

        const IOobject& io = *objects[fieldName];

        fields.set
        (
            nFields++,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobjectOption::MUST_READ,
                    IOobjectOption::AUTO_WRITE,
                    io.registerObject()
                ),
                mesh,
                readOldTime
            )
        );
    }

    if (nFields) Info<< endl;

    return fieldNames;
}


template<class GeoField, class Mesh>
Foam::wordList Foam::ReadFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar
)
{
    // Names of GeoField objects, sorted order.
    const wordList fieldNames(objects.sortedNames<GeoField>(syncPar));

    // Construct the fields - reading in consistent (master) order.
    fields.resize_null(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Reading " << GeoField::typeName << ':';
        }
        Info<< ' ' << fieldName;

        const IOobject& io = *objects[fieldName];

        fields.set
        (
            nFields++,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobjectOption::MUST_READ,
                    IOobjectOption::AUTO_WRITE,
                    io.registerObject()
                ),
                mesh
            )
        );
    }

    if (nFields) Info<< endl;

    return fieldNames;
}


template<class GeoField>
Foam::wordList Foam::ReadFields
(
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar
)
{
    // Names of GeoField objects, sorted order.
    const wordList fieldNames(objects.sortedNames<GeoField>(syncPar));

    // Construct the fields - reading in consistent (master) order.
    fields.resize_null(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Reading " << GeoField::typeName << ':';
        }
        Info<< ' ' << fieldName;

        const IOobject& io = *objects[fieldName];

        fields.set
        (
            nFields++,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobjectOption::MUST_READ,
                    IOobjectOption::AUTO_WRITE,
                    io.registerObject()
                )
            )
        );
    }

    if (nFields) Info<< endl;

    return fieldNames;
}


template<class GeoField>
void Foam::ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    objectRegistry& fieldsCache
)
{
    // Unload times that are no longer used
    {
        wordHashSet unusedTimes(fieldsCache.toc());
        unusedTimes.erase(timeNames);

        //Info<< "Unloading times " << unusedTimes << endl;

        for (const word& timeName : unusedTimes)
        {
            objectRegistry& timeCache =
                fieldsCache.lookupObjectRef<objectRegistry>(timeName);

            fieldsCache.checkOut(timeCache);
        }
    }


    // Load any new fields
    for (const word& timeName : timeNames)
    {
        // Create if not found
        if (!fieldsCache.found(timeName))
        {
            //Info<< "Creating registry for time " << timeName << endl;

            // Create objectRegistry if not found
            objectRegistry* timeCachePtr = new objectRegistry
            (
                IOobject
                (
                    timeName,
                    timeName,
                    fieldsCache,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::REGISTER
                )
            );
            timeCachePtr->store();
        }

        // Obtain cache for current time
        const objectRegistry& timeCache =
            fieldsCache.lookupObject<objectRegistry>(timeName);

        // Store field if not found
        if (!timeCache.found(fieldName))
        {
            //Info<< "Loading field " << fieldName
            //    << " for time " << timeName << endl;

            GeoField loadedFld
            (
                IOobject
                (
                    fieldName,
                    timeName,
                    mesh.thisDb(),
                    IOobjectOption::MUST_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::NO_REGISTER
                ),
                mesh
            );

            // Transfer to timeCache (new objectRegistry and store flag)
            GeoField* fldPtr = new GeoField
            (
                IOobject
                (
                    fieldName,
                    timeName,
                    timeCache,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::REGISTER
                ),
                loadedFld
            );
            fldPtr->store();
        }
    }
}


template<class GeoField>
void Foam::ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    const word& registryName
)
{
    ReadFields<GeoField>
    (
        fieldName,
        mesh,
        timeNames,
        const_cast<objectRegistry&>
        (
            mesh.thisDb().subRegistry(registryName, true)
        )
    );
}


template<class GeoFieldType, class NameMatchPredicate>
void Foam::readFields
(
    const typename GeoFieldType::Mesh& mesh,
    const IOobjectList& objects,
    const NameMatchPredicate& selectedFields,
    DynamicList<regIOobject*>& storedObjects
)
{
    // GeoField objects, sorted order. Not synchronised.
    const UPtrList<const IOobject> fieldObjects
    (
        objects.csorted<GeoFieldType>(selectedFields)
    );


    // pre-extend reserve
    storedObjects.reserve(storedObjects.size() + fieldObjects.size());

    label nFields = 0;

    for (const IOobject& io : fieldObjects)
    {
        if (!nFields)
        {
            Info<< "    " << GeoFieldType::typeName << ':';
        }
        Info<< ' ' << io.name();

        GeoFieldType* fieldPtr = new GeoFieldType
        (
            IOobject
            (
                io.name(),
                io.instance(),
                io.local(),
                io.db(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            mesh
        );
        fieldPtr->store();
        storedObjects.push_back(fieldPtr);

        ++nFields;
    }

    if (nFields) Info<< endl;
}


template<class UniformFieldType, class NameMatchPredicate>
void Foam::readUniformFields
(
    const IOobjectList& objects,
    const NameMatchPredicate& selectedFields,
    DynamicList<regIOobject*>& storedObjects,
    const bool syncPar
)
{
    // UniformField objects, sorted order, synchronised.
    const UPtrList<const IOobject> fieldObjects
    (
        objects.csorted<UniformFieldType>(selectedFields, syncPar)
    );

    // pre-extend reserve
    storedObjects.reserve(storedObjects.size() + fieldObjects.size());

    label nFields = 0;

    for (const IOobject& io : fieldObjects)
    {
        if (!nFields)
        {
            Info<< "    " << UniformFieldType::typeName << ':';
        }
        Info<< ' ' << io.name();

        UniformFieldType* fieldPtr = new UniformFieldType
        (
            IOobject
            (
                io.name(),
                io.instance(),
                io.local(),
                io.db(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            )
        );
        fieldPtr->store();
        storedObjects.push_back(fieldPtr);

        ++nFields;
    }

    if (nFields) Info<< endl;
}


template<class GeoFieldType, class NameMatchPredicate>
void Foam::readFields
(
    const typename GeoFieldType::Mesh& mesh,
    const IOobjectList& objects,
    const NameMatchPredicate& selectedFields,
    LIFOStack<regIOobject*>& storedObjects
)
{
    DynamicList<regIOobject*> newObjects;

    readFields<GeoFieldType, NameMatchPredicate>
    (
        mesh,
        objects,
        selectedFields,
        newObjects
    );

    // Transcribe from list to stack
    for (regIOobject* fieldPtr : newObjects)
    {
        storedObjects.push(fieldPtr);
    }
}


template<class UniformFieldType, class NameMatchPredicate>
void Foam::readUniformFields
(
    const IOobjectList& objects,
    const NameMatchPredicate& selectedFields,
    LIFOStack<regIOobject*>& storedObjects,
    const bool syncPar
)
{
    DynamicList<regIOobject*> newObjects;

    readUniformFields<UniformFieldType, NameMatchPredicate>
    (
        objects,
        selectedFields,
        newObjects,
        syncPar
    );

    // Transcribe from list to stack
    for (regIOobject* fieldPtr : newObjects)
    {
        storedObjects.push(fieldPtr);
    }
}


// ************************************************************************* //
