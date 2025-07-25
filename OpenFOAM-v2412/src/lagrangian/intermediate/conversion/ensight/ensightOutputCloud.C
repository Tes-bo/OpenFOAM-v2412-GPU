/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "ensightOutputCloud.H"
#include "fvMesh.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Binary output
static inline void writeMeasured_binary
(
    ensightFile& os,
    const UList<floatVector>& points
)
{
    for (const auto& p : points)
    {
        os.write(p.x());
        os.write(p.y());
        os.write(p.z());
    }
}

//- ASCII output. Id + position together
static inline label writeMeasured_ascii
(
    ensightFile& os,
    label pointId,
    const UList<floatVector>& points
)
{
    for (const auto& p : points)
    {
        os.writeInt(++pointId, 8);  // 1-index and an unusual width
        os.write(p.x());
        os.write(p.y());
        os.write(p.z());
        os.newline();
    }

    return pointId;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::ensightOutput::writeCloudPositions
(
    ensightFile& os,
    DynamicList<floatVector>& positions,
    const globalIndex& procAddr
)
{
    // Total number of parcels across all ranks
    const label nTotParcels = procAddr.totalSize();

    bool noCloud(!procAddr.totalSize());
    Pstream::broadcast(noCloud);

    if (UPstream::master())
    {
        os.beginParticleCoordinates(nTotParcels);
    }

    if (noCloud)
    {
        return false;  // All empty
    }

    if (UPstream::master())
    {
        const bool isBinaryOutput = (os.format() == IOstreamOption::BINARY);

        label parcelId = 0;

        if (isBinaryOutput)
        {
            // NB: binary write is Ensight6 - first ids, then positions

            // 1-index
            for (label id = 1; id <= nTotParcels; ++id)
            {
                os.write(id);
            }

            // Write master data
            writeMeasured_binary(os, positions);
        }
        else
        {
            // NB: ascii write is (id + position) together

            // Write master data
            parcelId = writeMeasured_ascii(os, parcelId, positions);
        }


        positions.clear();
        positions.reserve_nocopy(procAddr.maxNonLocalSize());

        // Receive and write
        for (const label proci : procAddr.subProcs())
        {
            const label procSize = procAddr.localSize(proci);

            if (procSize)
            {
                positions.resize_nocopy(procSize);

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    positions.data_bytes(),
                    positions.size_bytes()
                );

                if (isBinaryOutput)
                {
                    writeMeasured_binary(os, positions);
                }
                else
                {
                    parcelId = writeMeasured_ascii(os, parcelId, positions);
                }
            }
        }
    }
    else if (UPstream::is_subrank())
    {
        if (positions.size())
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                positions.cdata_bytes(),
                positions.size_bytes()
            );
        }
    }

    return true;
}


bool Foam::ensightOutput::writeCloudPositions
(
    ensightFile& os,
    DynamicList<floatVector>& positions
)
{
    return ensightOutput::writeCloudPositions
    (
        os,
        positions,
        // Gather sizes (offsets irrelevant)
        globalIndex(globalIndex::gatherOnly{}, positions.size())
    );
}


bool Foam::ensightOutput::writeCloudPositions
(
    ensightFile& os,
    const fvMesh& mesh,
    const word& cloudName,
    bool exists
)
{
    autoPtr<Cloud<passiveParticle>> parcelsPtr;

    if (exists)
    {
        parcelsPtr.reset(new Cloud<passiveParticle>(mesh, cloudName, false));
    }

    const label nLocalParcels
    (
        parcelsPtr ? parcelsPtr->size() : 0
    );

    // Gather sizes (offsets irrelevant)
    // and total number of parcels (all processes)
    const globalIndex procAddr(globalIndex::gatherOnly{}, nLocalParcels);

    // Extract positions from parcel.
    // Store as floatVector, since that is what Ensight will write anyhow

    DynamicList<floatVector> positions;
    positions.reserve(UPstream::master() ? procAddr.maxSize() : nLocalParcels);

    if (parcelsPtr)
    {
        const auto& parcels = *parcelsPtr;

        positions.resize_nocopy(parcels.size());  // same as nLocalParcels

        auto iter = positions.begin();

        if (std::is_same<float, vector::cmptType>::value)
        {
            for (const auto& p : parcels)
            {
                *iter = p.position();
                ++iter;
            }
        }
        else
        {
            for (const auto& p : parcels)
            {
                const vector pos(p.position());

                (*iter).x() = narrowFloat(pos.x());
                (*iter).y() = narrowFloat(pos.y());
                (*iter).z() = narrowFloat(pos.z());
                ++iter;
            }
        }

        parcelsPtr.reset(nullptr);
    }

    return ensightOutput::writeCloudPositions(os, positions, procAddr);
}


// ************************************************************************* //
