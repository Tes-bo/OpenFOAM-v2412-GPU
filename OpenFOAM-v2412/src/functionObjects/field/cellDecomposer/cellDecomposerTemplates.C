/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "fvMesh.H"
#include "emptyFvPatchFields.H"
#include "directFvPatchFieldMapper.H"
#include "mapPolyMesh.H"
//#include "polyPatch.H"
//#include "lduSchedule.H"
//#include "meshToMesh.H"

//template<class Type>
//void Foam::functionObjects::cellDecomposer::evaluateConstraintTypes
//(
//    GeometricField<Type, fvPatchField, volMesh>& fld
//) const
//{
//    auto& fldBf = fld.boundaryFieldRef();
//
//    const UPstream::commsTypes commsType = UPstream::defaultCommsType;
//    const label startOfRequests = UPstream::nRequests();
//
//    if
//    (
//        commsType == UPstream::commsTypes::blocking
//     || commsType == UPstream::commsTypes::nonBlocking
//    )
//    {
//        forAll(fldBf, patchi)
//        {
//            fvPatchField<Type>& tgtField = fldBf[patchi];
//
//            if
//            (
//                tgtField.type() == tgtField.patch().patch().type()
//             && polyPatch::constraintType(tgtField.patch().patch().type())
//            )
//            {
//                tgtField.initEvaluate(commsType);
//            }
//        }
//
//        // Wait for outstanding requests
//        if (commsType == UPstream::commsTypes::nonBlocking)
//        {
//            UPstream::waitRequests(startOfRequests);
//        }
//
//        forAll(fldBf, patchi)
//        {
//            fvPatchField<Type>& tgtField = fldBf[patchi];
//
//            if
//            (
//                tgtField.type() == tgtField.patch().patch().type()
//             && polyPatch::constraintType(tgtField.patch().patch().type())
//            )
//            {
//                tgtField.evaluate(commsType);
//            }
//        }
//    }
//    else if (commsType == UPstream::commsTypes::scheduled)
//    {
//        const lduSchedule& patchSchedule =
//            fld.mesh().globalData().patchSchedule();
//
//        for (const auto& schedEval : patchSchedule)
//        {
//            const label patchi = schedEval.patch;
//
//            fvPatchField<Type>& tgtField = fldBf[patchi];
//
//            if
//            (
//                tgtField.type() == tgtField.patch().patch().type()
//             && polyPatch::constraintType(tgtField.patch().patch().type())
//            )
//            {
//                if (schedEval.init)
//                {
//                    tgtField.initEvaluate(commsType);
//                }
//                else
//                {
//                    tgtField.evaluate(commsType);
//                }
//            }
//        }
//    }
//}

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
>
Foam::functionObjects::cellDecomposer::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const fvMesh& sMesh,
    const labelUList& patchMap,
    const labelUList& cellMap,
    const labelUList& faceMap,
    const bool allowUnmapped
) const
{
    // 1. Create the complete field with dummy patch fields
    PtrList<fvPatchField<Type>> patchFields(patchMap.size());

    forAll(patchFields, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces. Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] == -1)
        {
            patchFields.set
            (
                patchi,
                new emptyFvPatchField<Type>
                (
                    sMesh.boundary()[patchi],
                    fvPatchField<Type>::Internal::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    fvPatchFieldBase::calculatedType(),
                    sMesh.boundary()[patchi],
                    fvPatchField<Type>::Internal::null()
                )
            );
        }
    }

    auto tresult = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "subset"+vf.name(),
            sMesh.time().timeName(),
            sMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sMesh,
        vf.dimensions(),
        Field<Type>(vf.primitiveField(), cellMap),
        patchFields
    );
    auto& result = tresult.ref();
    result.oriented() = vf.oriented();


    // 2. Change the fvPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(bf, patchi)
    {
        const label basePatchId = patchMap[patchi];

        if (basePatchId != -1)
        {
            // Construct addressing
            const fvPatch& subPatch = sMesh.boundary()[patchi];
            const fvPatch& basePatch = vf.mesh().boundary()[basePatchId];
            const label baseStart = basePatch.start();
            const label baseSize = basePatch.size();

            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                const label baseFacei = faceMap[subPatch.start()+i];

                if (baseFacei >= baseStart && baseFacei < baseStart+baseSize)
                {
                    directAddressing[i] = baseFacei-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Leave up to
                    // fvPatchField
                    directAddressing[i] = -1;
                }
            }


            directFvPatchFieldMapper mapper(directAddressing);

            // allowUnmapped : special mode for if we do not want to be
            // warned for unmapped faces (e.g. from fvMeshDistribute).

            const bool hasUnmapped = mapper.hasUnmapped();
            if (allowUnmapped)
            {
                mapper.hasUnmapped() = false;
            }

            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[basePatchId],
                    subPatch,
                    result.internalField(),
                    mapper
                )
            );

            if (allowUnmapped && hasUnmapped)
            {
                // Set unmapped values to zeroGradient. This is the default
                // action for unmapped fvPatchFields. Note that this bypasses
                // any special logic for handling unmapped fvPatchFields but
                // since this is only used inside fvMeshDistribute ...

                tmp<Field<Type>> tfld(bf[patchi].patchInternalField());
                const Field<Type>& fld = tfld();

                Field<Type> value(bf[patchi]);
                forAll(directAddressing, i)
                {
                    if (directAddressing[i] == -1)
                    {
                        value[i] = fld[i];
                    }
                }
                bf[patchi].fvPatchField<Type>::operator=(value);
            }
        }
    }

    return tresult;
}


template<class Type>
bool Foam::functionObjects::cellDecomposer::mapFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const fvMesh& mapRegion =
        this->mesh_.time().lookupObject<fvMesh>(mapRegion_);

    const labelList patchMap(identity(mapRegion.boundaryMesh().size()));

    const wordList fieldNames
    (
        this->mesh_.sortedNames<VolFieldType>(fieldNames_)
    );

    const bool processed = !fieldNames.empty();

    for (const word& fieldName : fieldNames)
    {
        const VolFieldType& field = lookupObject<VolFieldType>(fieldName);

        auto* mapFieldPtr = mapRegion.getObjectPtr<VolFieldType>(fieldName);

        if (!mapFieldPtr)
        {
            mapFieldPtr = new VolFieldType
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    mapRegion,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                mapRegion,
                dimensioned<Type>(field.dimensions(), Zero)
            );

            mapFieldPtr->store();
        }

        auto& mappedField = *mapFieldPtr;

        mappedField = interpolate
        (
            field,
            mapRegion,
            patchMap,
            mapPtr_().cellMap(),
            mapPtr_().faceMap(),
            false             //allowUnmapped
        );
        Log << "    " << fieldName << ": interpolated";

        //evaluateConstraintTypes(mappedField);
    }

    return processed;
}


template<class Type>
bool Foam::functionObjects::cellDecomposer::writeFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const fvMesh& mapRegion =
        this->mesh_.time().lookupObject<fvMesh>(mapRegion_);

    const wordList fieldNames
    (
        this->mesh_.sortedNames<VolFieldType>(fieldNames_)
    );

    const bool processed = !fieldNames.empty();

    for (const word& fieldName : fieldNames)
    {
        const VolFieldType& mappedField =
            mapRegion.template lookupObject<VolFieldType>(fieldName);

        mappedField.write();

        Log << "    " << fieldName << ": written";
    }

    return processed;
}


// ************************************************************************* //
