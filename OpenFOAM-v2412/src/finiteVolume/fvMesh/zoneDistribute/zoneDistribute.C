/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "zoneDistribute.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneDistribute, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneDistribute::zoneDistribute(const fvMesh& mesh)
:
    MeshObject_type(mesh),
    stencil_(zoneCPCStencil::New(mesh)),
    globalNumbering_(stencil_.globalNumbering()),
    send_(UPstream::nProcs()),
    pBufs_(UPstream::commsTypes::nonBlocking)
{
    // Don't clear storage on persistent buffer
    pBufs_.allowClearRecv(false);
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::zoneDistribute& Foam::zoneDistribute::New(const fvMesh& mesh)
{
    auto* ptr = mesh.thisDb().getObjectPtr<zoneDistribute>("zoneDistribute");

    if (!ptr)
    {
        ptr = new zoneDistribute(mesh);
        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneDistribute::updateStencil(const boolList& zone)
{
    zoneCPCStencil::New(mesh_).updateStencil(zone);
}


void Foam::zoneDistribute::setUpCommforZone
(
    const boolList& zone,
    bool updateStencil
)
{
    zoneCPCStencil& stencil = zoneCPCStencil::New(mesh_);

    if (updateStencil)
    {
        stencil.updateStencil(zone);
    }

    if (UPstream::parRun())
    {
        List<labelHashSet> needed(UPstream::nProcs());

        // Bin according to originating (sending) processor
        for (const label celli : stencil.needsComm())
        {
            if (zone[celli])
            {
                for (const label gblIdx : stencil_[celli])
                {
                    const label proci = globalNumbering_.whichProcID(gblIdx);

                    if (proci != Pstream::myProcNo())
                    {
                        needed[proci].insert(gblIdx);
                    }
                }
            }
        }

        // Stream the send data into PstreamBuffers,
        // which we also use to track the current topology.

        pBufs_.clear();

        for (const int proci : pBufs_.allProcs())
        {
            const auto& indices = needed[proci];

            if (proci != UPstream::myProcNo() && !indices.empty())
            {
                // Serialize as List
                UOPstream toProc(proci, pBufs_);
                toProc << indices.sortedToc();
            }
        }

        pBufs_.finishedSends(sendConnections_, sendProcs_, recvProcs_);

        for (const int proci : pBufs_.allProcs())
        {
            send_[proci].clear();

            if (proci != UPstream::myProcNo() && pBufs_.recvDataCount(proci))
            {
                UIPstream fromProc(proci, pBufs_);
                fromProc >> send_[proci];
            }
        }
    }
}


// ************************************************************************* //
