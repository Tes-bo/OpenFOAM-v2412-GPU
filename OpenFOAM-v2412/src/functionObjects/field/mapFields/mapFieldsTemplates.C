/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
#include "polyPatch.H"
#include "lduSchedule.H"
#include "meshToMesh.H"

template<class Type>
void Foam::functionObjects::mapFields::evaluateConstraintTypes
(
    GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    auto& bfld = fld.boundaryFieldRef();

    const UPstream::commsTypes commsType = UPstream::defaultCommsType;

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : bfld)
        {
            if
            (
                pfld.type() == pfld.patch().patch().type()
             && polyPatch::constraintType(pfld.patch().patch().type())
            )
            {
                pfld.initEvaluate(commsType);
            }
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : bfld)
        {
            if
            (
                pfld.type() == pfld.patch().patch().type()
             && polyPatch::constraintType(pfld.patch().patch().type())
            )
            {
                pfld.evaluate(commsType);
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = bfld[patchi];

            if
            (
                pfld.type() == pfld.patch().patch().type()
             && polyPatch::constraintType(pfld.patch().patch().type())
            )
            {
                if (schedEval.init)
                {
                    pfld.initEvaluate(commsType);
                }
                else
                {
                    pfld.evaluate(commsType);
                }
            }
        }
    }
}


template<class Type>
bool Foam::functionObjects::mapFields::mapFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const fvMesh& mapRegion = mapRegionPtr_();

    wordList fieldNames(this->mesh_.sortedNames<VolFieldType>(fieldNames_));

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

        mappedField = interpPtr_->mapTgtToSrc(field);

        Log << "    " << fieldName << ": interpolated\n";

        evaluateConstraintTypes(mappedField);
    }

    return processed;
}


template<class Type>
bool Foam::functionObjects::mapFields::writeFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const fvMesh& mapRegion = mapRegionPtr_();

    wordList fieldNames(this->mesh_.sortedNames<VolFieldType>(fieldNames_));

    const bool processed = !fieldNames.empty();

    for (const word& fieldName : fieldNames)
    {
        const VolFieldType& mappedField =
            mapRegion.template lookupObject<VolFieldType>(fieldName);

        mappedField.write();

        Log << "    " << fieldName << ": written\n";
    }

    return processed;
}


// ************************************************************************* //
