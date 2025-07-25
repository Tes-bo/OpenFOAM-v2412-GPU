/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void Foam::fieldsDistributor::readField
(
    const IOobject& io,
    const typename GeoField::Mesh& mesh,
    const label idx,
    PtrList<GeoField>& fields
)
{
    fields.emplace_set(idx, io, mesh);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fieldsDistributor::readField
(
    const IOobject& io,
    const typename GeoMesh::Mesh& mesh,
    const label idx,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields
)
{
    fields.emplace_set(idx, io, mesh, false);  // readOldTime = false
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fieldsDistributor::readFields
(
    const typename GeoMesh::Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
    const bool readOldTime
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    // GeoField fields - sorted for consistent order on all processors
    const UPtrList<const IOobject> fieldObjects(objects.csorted<GeoField>());

    // Construct the fields
    fields.resize_null(fieldObjects.size());

    forAll(fieldObjects, i)
    {
        fields.emplace_set(i, fieldObjects[i], mesh, readOldTime);
    }
}


template<class Mesh, class GeoField>
void Foam::fieldsDistributor::readFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields
)
{
    // GeoField fields - sorted for consistent order on all processors
    const UPtrList<const IOobject> fieldObjects(objects.csorted<GeoField>());

    // Construct the fields
    fields.resize_null(fieldObjects.size());

    forAll(fieldObjects, i)
    {
        fields.emplace_set(i, fieldObjects[i], mesh);
    }
}


template<class BoolListType, class GeoField, class MeshSubsetter>
void Foam::fieldsDistributor::readFieldsImpl
(
    refPtr<fileOperation>* readHandlerPtr,  // Can be nullptr
    const BoolListType& haveMeshOnProc,
    const MeshSubsetter* subsetter,
    const typename GeoField::Mesh& mesh,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields,
    const bool deregister
)
{
    // Get my objects of type
    IOobjectList objects(allObjects.lookupClass<GeoField>());

    // Check that we all have all objects
    wordList objectNames = objects.sortedNames();

    // Get master names
    wordList masterNames(objectNames);
    Pstream::broadcast(masterNames);

    if
    (
        haveMeshOnProc.test(UPstream::myProcNo())
     && objectNames != masterNames
    )
    {
        FatalErrorInFunction
            << "Objects not synchronised across processors." << nl
            << "Master has " << flatOutput(masterNames) << nl
            << "Processor " << UPstream::myProcNo()
            << " has " << flatOutput(objectNames)
            << exit(FatalError);
    }

    #ifdef FULLDEBUG
    {
        // A bit more checking - this may not be strictly correct.
        // The master must know about everyone, but the others really
        // only need to know about themselves.
        // Can broadcast decompose = yes/no from master

        bitSet localValues(haveMeshOnProc);
        bitSet masterValues(localValues);
        Pstream::broadcast(masterValues);

        localValues ^= masterValues;

        if (localValues.any())
        {
            FatalErrorInFunction
                << "haveMeshOnProc not synchronised across processors." << nl
                << "Processor " << UPstream::myProcNo()
                << " differs at these positions: "
                << flatOutput(localValues.sortedToc()) << nl
                << exit(FatalError);
         }
    }
    #endif


    fields.resize_null(masterNames.size());

    if (fields.empty())
    {
        if (deregister)
        {
            // Extra safety - remove all such types
            for
            (
                const GeoField& fld
              : mesh.thisDb().objectRegistry::template csorted<GeoField>()
            )
            {
                if (!fld.ownedByRegistry())
                {
                    const_cast<GeoField&>(fld).checkOut();
                }
            }
        }

        // Early exit
        return;
    }


    // We are decomposing if none of the subprocs has a mesh
    bool decompose = true;
    if (UPstream::master())
    {
        for (const int proci : UPstream::subProcs())
        {
            if (haveMeshOnProc.test(proci))
            {
                decompose = false;
                break;
            }
        }
    }
    Pstream::broadcast(decompose);


    if (decompose && UPstream::master())
    {
        const bool oldParRun = UPstream::parRun(false);

        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt(IOobject::AUTO_WRITE);

            // Load field (but not oldTime)
            readField(io, mesh, i, fields);
        }

        UPstream::parRun(oldParRun);
    }
    else if
    (
        !decompose
     &&
        // Has read-handler : use it to decide if reading is possible
        // No  read-handler : decide based on the presence of a mesh
        (
            readHandlerPtr
          ? readHandlerPtr->good()
          : haveMeshOnProc.test(UPstream::myProcNo())
        )
    )
    {
        const label oldWorldComm = UPstream::worldComm;
        refPtr<fileOperation> oldHandler;

        if (readHandlerPtr)
        {
            // Swap read fileHandler for read fields
            oldHandler = fileOperation::fileHandler(*readHandlerPtr);
            UPstream::commWorld(fileHandler().comm());
        }

        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt(IOobject::AUTO_WRITE);

            // Load field (but not oldTime)
            readField(io, mesh, i, fields);
        }

        if (readHandlerPtr)
        {
            // Restore fileHandler
            *readHandlerPtr = fileOperation::fileHandler(oldHandler);
            UPstream::commWorld(oldWorldComm);
        }
    }


    // Missing fields on any processors?
    // - construct from dictionary

    PtrList<dictionary> fieldDicts;

    if (UPstream::master())
    {
        // Broadcast zero sized fields everywhere (if needed)
        // Send like a list of dictionaries

        OPBstream toProcs(UPstream::worldComm);

        const label nDicts = (subsetter ? fields.size() : label(0));

        toProcs << nDicts << token::BEGIN_LIST;  // Begin list

        if (nDicts && subsetter)
        {
            // Disable communication for interpolate() method
            const bool oldParRun = UPstream::parRun(false);

            for (const auto& fld : fields)
            {
                tmp<GeoField> tsubfld = subsetter->interpolate(fld);

                // Surround each with {} as dictionary entry
                toProcs.beginBlock();
                toProcs << tsubfld();
                toProcs.endBlock();
            }

            UPstream::parRun(oldParRun);  // Restore state
        }

        toProcs << token::END_LIST << token::NL;  // End list
    }
    else
    {
        // Receive the broadcast...
        IPBstream fromMaster(UPstream::worldComm);

        // But only consume where needed...
        if (!haveMeshOnProc.test(UPstream::myProcNo()))
        {
            fromMaster >> fieldDicts;
        }
    }


    // Use the received dictionaries (if any) to create missing fields.

    // Disable communication when constructing from dictionary
    const bool oldParRun = UPstream::parRun(false);

    IOobject noreadIO
    (
        "none",
        mesh.time().timeName(),
        mesh.thisDb(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    );

    forAll(fieldDicts, i)
    {
        noreadIO.resetHeader(masterNames[i]);

        fields.emplace_set(i, noreadIO, mesh, fieldDicts[i]);
    }

    UPstream::parRun(oldParRun);  // Restore any changes


    // Finally. Can checkOut of registry as required
    if (deregister)
    {
        /// Info<< "De-registering fields:";
        for (auto& fld : fields)
        {
            /// Info<< "  " << fld.name();

            // Ensure it is not destroyed by polyMesh deletion
            fld.checkOut();
        }
        /// Info<< nl;

        // Extra safety - remove all such types
        for
        (
            const GeoField& fld
          : mesh.thisDb().objectRegistry::template csorted<GeoField>()
        )
        {
            if (!fld.ownedByRegistry())
            {
                const_cast<GeoField&>(fld).checkOut();
            }
        }
    }
}


template<class GeoField, class MeshSubsetter>
void Foam::fieldsDistributor::readFields
(
    const bitSet& haveMeshOnProc,
    const MeshSubsetter* subsetter,
    const typename GeoField::Mesh& mesh,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields,
    const bool deregister
)
{
    readFieldsImpl
    (
        nullptr,  // readHandler
        haveMeshOnProc,
        subsetter,

        mesh,
        allObjects,
        fields,
        deregister
    );
}


template<class GeoField, class MeshSubsetter>
void Foam::fieldsDistributor::readFields
(
    const boolUList& haveMeshOnProc,
    const MeshSubsetter* subsetter,
    const typename GeoField::Mesh& mesh,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields,
    const bool deregister
)
{
    readFieldsImpl
    (
        nullptr,  // readHandler
        haveMeshOnProc,
        subsetter,

        mesh,
        allObjects,
        fields,
        deregister
    );
}


template<class GeoField, class MeshSubsetter>
void Foam::fieldsDistributor::readFields
(
    const boolUList& haveMeshOnProc,
    const typename GeoField::Mesh& mesh,
    const autoPtr<MeshSubsetter>& subsetter,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields,
    const bool deregister
)
{
    readFieldsImpl
    (
        nullptr,  // readHandler
        haveMeshOnProc,
        subsetter.get(),

        mesh,
        allObjects,
        fields,
        deregister
    );
}


template<class GeoField, class MeshSubsetter>
void Foam::fieldsDistributor::readFields
(
    const boolUList& haveMeshOnProc,
    refPtr<fileOperation>& readHandler,
    const typename GeoField::Mesh& mesh,
    const autoPtr<MeshSubsetter>& subsetter,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields,
    const bool deregister
)
{
    readFieldsImpl
    (
       &readHandler,
        haveMeshOnProc,
        subsetter.get(),

        mesh,
        allObjects,
        fields,
        deregister
    );
}


// ************************************************************************* //
