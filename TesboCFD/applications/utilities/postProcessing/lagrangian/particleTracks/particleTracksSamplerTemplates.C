/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "objectRegistry.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::particleTracksSampler::createTrackField
(
    const UList<Type>& values,
    List<DynamicList<Type>>& trackValues
) const
{
    List<Type> allValues(cloudGather_.gather(values));

    const label nTracks = trackValues.size();

    forAll(allValues, i)
    {
        const label globalId =
            origParcelAddr_.toGlobal(origProcIds_[i], origParcelIds_[i]);

        if (globalId % stride_ == 0)
        {
            const label trackId = globalId/stride_;

            if
            (
                trackId < nTracks
             && trackValues[trackId].size() < maxPositions_
            )
            {
                trackValues[trackId].append(allValues[i]);
            }
        }
    }
}


template<class Type>
Foam::label Foam::particleTracksSampler::setTrackFields
(
    const objectRegistry& obr,
    HashTable<List<DynamicList<Type>>>& fieldTable
) const
{
    // First time - obtain field names and populate storage locations
    if (fieldTable.empty())
    {
        wordList fieldNames = obr.names<IOField<Type>>();

        Pstream::combineReduce(fieldNames, ListOps::uniqueEqOp<word>());

        for (const word& fieldName : fieldNames)
        {
            fieldTable(fieldName).resize(nTracks());
        }
    }

    // Process in parallel-consistent order
    for (const word& fieldName : fieldTable.sortedToc())
    {
        auto& output = fieldTable[fieldName];

        const auto* ptr = obr.cfindObject<IOField<Type>>(fieldName);

        this->createTrackField
        (
            ptr ? static_cast<const UList<Type>&>(*ptr) : UList<Type>::null(),
            output
        );
    }

    return fieldTable.size();
}


// ************************************************************************* //
