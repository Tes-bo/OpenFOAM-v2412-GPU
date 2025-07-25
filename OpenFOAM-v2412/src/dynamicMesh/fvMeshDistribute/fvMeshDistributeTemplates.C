/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType, class ZoneMesh>
void Foam::fvMeshDistribute::reorderZones
(
    const wordList& zoneNames,
    ZoneMesh& zones
)
{
    zones.clearAddressing();

    // Shift old ones to new position
    UPtrList<ZoneType> newZonePtrs(zoneNames.size());
    forAll(zones, zonei)
    {
        auto* zonePtr = zones.get(zonei);
        if (!zonePtr)
        {
            FatalErrorInFunction << "Problem with zones " << zones.names()
                << exit(FatalError);
        }
        const label newIndex = zoneNames.find(zonePtr->name());
        zonePtr->index() = newIndex;
        newZonePtrs.set(newIndex, zonePtr);
    }

    // Add empty zones for unknown ones
    forAll(newZonePtrs, i)
    {
        if (!newZonePtrs.get(i))
        {
            newZonePtrs.set
            (
                i,
                new ZoneType
                (
                    zoneNames[i],
                    i,
                    zones
                )
            );
        }
    }

    // Transfer
    zones.swap(newZonePtrs);
}


template<class GeoField>
void Foam::fvMeshDistribute::printIntFieldInfo(const fvMesh& mesh)
{
    typedef GeometricField
    <
        typename GeoField::value_type,
        fvPatchField,
        volMesh
    > excludeType;

    for (const GeoField& field : mesh.objectRegistry::csorted<GeoField>())
    {
        if (!isA<excludeType>(field))
        {
            Pout<< "Field:" << field.name() << " size:" << field.size()
                //<< " value:" << field
                << endl;
        }
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::printFieldInfo(const fvMesh& mesh)
{
    for (const GeoField& field : mesh.objectRegistry::csorted<GeoField>())
    {
        Pout<< "Field:" << field.name() << " size:" << field.size()
            //<< " value:" << field
            << endl;

        for (const auto& patchFld : field.boundaryField())
        {
            Pout<< "    " << patchFld.patch().index()
                << ' ' << patchFld.patch().name()
                << ' ' << patchFld.type()
                << ' ' << patchFld.size()
                << nl;
        }
    }
}


template<class T, class Mesh>
void Foam::fvMeshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvsPatchField, T>>& bflds
) const
{
    // Save whole boundary field

    typedef GeometricField<T, fvsPatchField, Mesh> fldType;

    const UPtrList<const fldType> flds
    (
        mesh_.objectRegistry::csorted<fldType>()
    );

    bflds.resize_null(flds.size());

    label fieldi = 0;
    for (const fldType& fld : flds)
    {
        bflds.set(fieldi, fld.boundaryField().clone());

        ++fieldi;
    }
}


template<class T, class Mesh>
void Foam::fvMeshDistribute::mapBoundaryFields
(
    const mapPolyMesh& map,
    const PtrList<FieldField<fvsPatchField, T>>& oldBflds
)
{
    // Map boundary field

    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvsPatchField, Mesh> fldType;

    UPtrList<fldType> flds
    (
        mesh_.objectRegistry::sorted<fldType>()
    );

    if (flds.size() != oldBflds.size())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    forAll(flds, fieldi)
    {
        auto& bfld = flds[fieldi].boundaryFieldRef();
        const auto& oldBfld = oldBflds[fieldi];

        // Pull from old boundary field into bfld.

        forAll(bfld, patchi)
        {
            fvsPatchField<T>& patchFld = bfld[patchi];
            label facei = patchFld.patch().start();

            forAll(patchFld, i)
            {
                label oldFacei = faceMap[facei++];

                // Find patch and local patch face oldFacei was in.
                forAll(oldPatchStarts, oldPatchi)
                {
                    label oldLocalI = oldFacei - oldPatchStarts[oldPatchi];

                    if (oldLocalI >= 0 && oldLocalI < oldBfld[oldPatchi].size())
                    {
                        patchFld[i] = oldBfld[oldPatchi][oldLocalI];
                    }
                }
            }
        }
    }
}


template<class T>
void Foam::fvMeshDistribute::saveInternalFields
(
    PtrList<Field<T>>& iflds
) const
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> fldType;

    const UPtrList<const fldType> fields
    (
        mesh_.objectRegistry::csorted<fldType>()
    );

    iflds.resize_null(fields.size());

    forAll(fields, fieldi)
    {
        iflds.set(fieldi, fields[fieldi].primitiveField().clone());
    }
}


template<class T>
void Foam::fvMeshDistribute::mapExposedFaces
(
    const mapPolyMesh& map,
    const PtrList<Field<T>>& oldFlds
)
{
    // Set boundary values of exposed internal faces

    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvsPatchField, surfaceMesh> fldType;

    UPtrList<fldType> flds
    (
        mesh_.objectRegistry::sorted<fldType>()
    );

    if (flds.size() != oldFlds.size())
    {
        FatalErrorInFunction
            << "problem"
            << abort(FatalError);
    }


    forAll(flds, fieldi)
    {
        auto& fld = flds[fieldi];
        const auto& oldInternal = oldFlds[fieldi];

        const bool oriented = fld.is_oriented();

        auto& bfld = fld.boundaryFieldRef();


        // Pull from old internal field into bfld.

        forAll(bfld, patchi)
        {
            fvsPatchField<T>& patchFld = bfld[patchi];

            forAll(patchFld, i)
            {
                const label faceI = patchFld.patch().start()+i;

                label oldFaceI = faceMap[faceI];

                if (oldFaceI < oldInternal.size())
                {
                    patchFld[i] = oldInternal[oldFaceI];

                    if (oriented && map.flipFaceFlux().found(faceI))
                    {
                        patchFld[i] = flipOp()(patchFld[i]);
                    }
                }
            }
        }
    }
}


template<class GeoField, class PatchFieldType>
void Foam::fvMeshDistribute::initPatchFields
(
    const typename GeoField::value_type& initVal
)
{
    // Init patch fields of certain type
    // - field order is irrelevant

    for (GeoField& fld : mesh_.objectRegistry::objects<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        forAll(bfld, patchi)
        {
            if (isA<PatchFieldType>(bfld[patchi]))
            {
                bfld[patchi] == initVal;
            }
        }
    }
}


//template<class GeoField>
//void Foam::fvMeshDistribute::correctBoundaryConditions()
//{
//    // CorrectBoundaryConditions patch fields of certain type
//
//    for (GeoField& fld : mesh_.objectRegistry::sorted<GeoField>())
//    {
//        fld.correctBoundaryConditions();
//    }
//}


template<class GeoField>
void Foam::fvMeshDistribute::getFieldNames
(
    const fvMesh& mesh,
    HashTable<wordList>& allFieldNames,
    const word& excludeType,
    const bool syncPar
)
{
    wordList& list = allFieldNames(GeoField::typeName);
    list = mesh.sortedNames<GeoField>();

    if (!excludeType.empty())
    {
        const wordList& excludeList =
            allFieldNames.lookup(excludeType, wordList::null());

        if (!excludeList.empty())
        {
            DynamicList<word> newList(list.size());
            for (const auto& name : list)
            {
                if (!excludeList.contains(name))
                {
                    newList.push_back(name);
                }
            }
            if (newList.size() < list.size())
            {
                list = std::move(newList);
            }
        }
    }


    // Check all procs have same names
    if (syncPar && Pstream::parRun())
    {
        // Check and report error(s) on master
        // - don't need indexing on master itself

        const globalIndex procAddr(globalIndex::gatherNonLocal{}, list.size());

        const wordList allNames(procAddr.gather(list));

        // Automatically restricted to master
        for (const int proci : procAddr.subProcs())
        {
            const auto procNames(allNames.slice(procAddr.range(proci)));

            if (procNames != list)
            {
                FatalErrorInFunction
                    << "When checking for equal " << GeoField::typeName
                    << " :" << nl
                    << "processor0 has:" << list << nl
                    << "processor" << proci << " has:" << procNames << nl
                    << GeoField::typeName
                    << " need to be synchronised on all processors."
                    << exit(FatalError);
                break;
            }
        }
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::sendFields
(
    const label domain,
    const HashTable<wordList>& allFieldNames,
    const fvMeshSubset& subsetter,
    Ostream& toNbr
)
{
    // Send fields. Note order supplied so we can receive in exactly the same
    // order.
    // Note that field gets written as entry in dictionary so we
    // can construct from subdictionary.
    // (since otherwise the reading as-a-dictionary mixes up entries from
    // consecutive fields)
    // The dictionary constructed is:
    //  volScalarField
    //  {
    //      p {internalField ..; boundaryField ..;}
    //      k {internalField ..; boundaryField ..;}
    //  }
    //  volVectorField
    //  {
    //      U {internalField ...  }
    //  }

    // volVectorField {U {internalField ..; boundaryField ..;}}

    const wordList& fieldNames =
        allFieldNames.lookup(GeoField::typeName, wordList::null());

    toNbr << GeoField::typeName << token::NL << token::BEGIN_BLOCK << token::NL;

    for (const word& fieldName : fieldNames)
    {
        if (debug)
        {
            Pout<< "Subsetting " << GeoField::typeName
                << " field " << fieldName
                << " for domain:" << domain << endl;
        }

        // Send all fieldNames. This has to be exactly the same set as is
        // being received!
        const GeoField& fld =
            subsetter.baseMesh().lookupObject<GeoField>(fieldName);

        // Note: use subsetter to get sub field. Override default behaviour
        //       to warn for unset fields since they will be reset later on
        tmp<GeoField> tsubfld = subsetter.interpolate(fld, true);

        toNbr
            << fieldName << token::NL << token::BEGIN_BLOCK
            << tsubfld
            << token::NL << token::END_BLOCK << token::NL;
    }
    toNbr << token::END_BLOCK << token::NL;
}


template<class GeoField>
void Foam::fvMeshDistribute::receiveFields
(
    const label domain,
    const HashTable<wordList>& allFieldNames,
    fvMesh& mesh,
    PtrList<GeoField>& fields,
    const dictionary& allFieldsDict
)
{
    // Opposite of sendFields

    const wordList& fieldNames =
        allFieldNames.lookup(GeoField::typeName, wordList::null());

    const dictionary& fieldDicts =
        allFieldsDict.subDict(GeoField::typeName);


    if (debug)
    {
        Pout<< "Receiving:" << GeoField::typeName
            << " fields:" << fieldNames
            << " from domain:" << domain << endl;
    }

    fields.resize(fieldNames.size());

    label fieldi = 0;
    for (const word& fieldName : fieldNames)
    {
        if (debug)
        {
            Pout<< "Constructing type:"  << GeoField::typeName
                << " field:" << fieldName
                << " from domain:" << domain << endl;
        }

        fields.set
        (
            fieldi++,
            new GeoField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                fieldDicts.subDict(fieldName)
            )
        );
    }
}


// ************************************************************************* //
