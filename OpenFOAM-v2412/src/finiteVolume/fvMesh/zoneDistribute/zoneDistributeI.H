/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "DynamicField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<typename Type>
Type Foam::zoneDistribute::getLocalValue
(
    const VolumeField<Type>& phi,
    const label localIdx
) const
{
    if (localIdx < mesh_.nCells()) // internal: cellI
    {
        return phi[localIdx];
    }

    return faceValue(phi,localIdx);
}


template<typename Type>
Type Foam::zoneDistribute::faceValue
(
    const VolumeField<Type>& phi,
    const label localIdx
) const
{
    const label faceI = localIdx + mesh_.nInternalFaces() - mesh_.nCells();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Boundary face. Find out which face of which patch
    const label patchI = pbm.whichPatch(faceI);

    if (patchI < 0 || patchI >= pbm.size())
    {
       FatalErrorInFunction
            << "Cannot find patch for face " << faceI
            << abort(FatalError);
    }

    const polyPatch& pp = pbm[patchI];

    const label patchFaceI = pp.whichFace(faceI);

    return phi.boundaryField()[patchI][patchFaceI];
}


template<typename Type>
Type Foam::zoneDistribute::getValue
(
    const VolumeField<Type>& phi,
    const Map<Type>& valuesFromOtherProc,
    const label gblIdx
) const
{
    if (globalNumbering_.isLocal(gblIdx))
    {
        const label localIdx = globalNumbering_.toLocal(gblIdx);
        return getLocalValue(phi,localIdx);
    }
    else
    {
        // From other proc
        return valuesFromOtherProc[gblIdx];
    }
}


template<typename Type>
Foam::Map<Foam::Field<Type>> Foam::zoneDistribute::getFields
(
    const boolList& zone,
    const VolumeField<Type>& phi
)
{
    if (zone.size() != phi.size())
    {
        FatalErrorInFunction
            << "size of zone: " << zone.size()
            << "size of phi:" <<  phi.size()
            << "do not match. Did the mesh change?"
            << exit(FatalError);
    }


    // Get values from other proc
    Map<Type> neiValues = getDatafromOtherProc(zone, phi);

    Map<Field<Type>> stencilWithValues;

    DynamicField<Type> tmpField(128);

    forAll(zone, celli)
    {
        if (zone[celli])
        {
            tmpField.clear();

            for (const label gblIdx : stencil_[celli])
            {
                tmpField.append(getValue(phi,neiValues,gblIdx));
            }

            stencilWithValues.emplace(celli, tmpField);
        }
    }

    return stencilWithValues;
}


template<typename Type>
Foam::Map<Type> Foam::zoneDistribute::getDatafromOtherProc
(
    const boolList& zone,
    const VolumeField<Type>& phi
)
{
    if (zone.size() != phi.size())
    {
        FatalErrorInFunction
            << "size of zone: " << zone.size()
            << "size of phi:" <<  phi.size()
            << "do not match. Did the mesh change?"
            << exit(FatalError);
    }


    // Get values from other proc
    Map<Type> neiValues;

    if (UPstream::parRun())
    {
        if (sendConnections_.empty())
        {
            WarningInFunction
                << "The send/recv connections not initialized - "
                << "likely that setUpCommforZone() was not called"
                << endl;
            // But don't exit/abort for now
        }

        // Stream the send data into PstreamBuffers,
        // which we also use to track the current topology.

        pBufs_.clear();

        for (const int proci : pBufs_.allProcs())
        {
            const auto& indices = send_[proci];

            if (proci != UPstream::myProcNo() && !indices.empty())
            {
                // Serialize as Map
                Map<Type> sendValues(2*indices.size());

                for (const label sendIdx : indices)
                {
                    sendValues.insert
                    (
                        sendIdx,
                        getLocalValue(phi, globalNumbering_.toLocal(sendIdx))
                    );
                }

                UOPstream toProc(proci, pBufs_);
                toProc << sendValues;
            }
        }

        pBufs_.finishedSends(sendConnections_, sendProcs_, recvProcs_);

        for (const int proci : pBufs_.allProcs())
        {
            if (proci != UPstream::myProcNo() && pBufs_.recvDataCount(proci))
            {
                UIPstream fromProc(proci, pBufs_);
                Map<Type> tmpValues(fromProc);

                neiValues += tmpValues;
            }
        }
    }

    return neiValues;
}


// ************************************************************************* //
