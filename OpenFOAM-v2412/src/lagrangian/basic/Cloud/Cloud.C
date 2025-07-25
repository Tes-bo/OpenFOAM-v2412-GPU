/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
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

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamBuffers.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::checkPatches() const
{
    for (const polyPatch& pp : polyMesh_.boundaryMesh())
    {
        const auto* camipp = isA<cyclicAMIPolyPatch>(pp);

        if (camipp && camipp->owner() && camipp->AMI().distributed())
        {
            FatalErrorInFunction
                << "Particle tracking across AMI patches is only currently "
                << "supported for cases where the AMI patches reside on a "
                << "single processor" << abort(FatalError);
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const Foam::zero,
    const word& cloudName
)
:
    cloud(pMesh, cloudName),
    polyMesh_(pMesh),
    geometryType_(cloud::geometryType::COORDINATES)
{
    checkPatches();

    (void)polyMesh_.tetBasePtIs();
    (void)polyMesh_.oldCellCentres();
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    Cloud<ParticleType>(pMesh, Foam::zero{}, cloudName)
{
    if (particles.size())
    {
        IDLList<ParticleType>::operator=(particles);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::addParticle(ParticleType* pPtr)
{
    this->append(pPtr);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteParticle(ParticleType& p)
{
    delete(this->remove(&p));
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteLostParticles()
{
    for (ParticleType& p : *this)
    {
        if (p.cell() == -1)
        {
            WarningInFunction
                << "deleting lost particle at position " << p.position()
                << endl;

            deleteParticle(p);
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::cloudReset(const Cloud<ParticleType>& c)
{
    // Reset particle count and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    ParticleType::particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::Cloud<ParticleType>::move
(
    TrackCloudType& cloud,
    typename ParticleType::trackingData& td,
    const scalar trackTime
)
{
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const globalMeshData& pData = polyMesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData.topology().procNeighbours();

    // Initialise the stepFraction moved for the particles
    for (ParticleType& p : *this)
    {
        p.reset();
    }

    // Clear the global positions as these are about to change
    globalPositionsPtr_.clear();


    // For v2112 and earlier: pre-assembled lists of particles
    // to be transferred and target patch on a per processor basis.
    // Apart from memory overhead of assembling the lists this adds
    // allocations/de-allocation when building linked-lists.

    // Now stream particle transfer tuples directly into PstreamBuffers.
    // Use a local cache of UOPstream wrappers for the formatters
    // (since there are potentially many particles being shifted about).


    // Allocate transfer buffers,
    // automatic clearStorage when UIPstream closes is disabled.
    PstreamBuffers pBufs;
    pBufs.allowClearRecv(false);

    // Cache of opened UOPstream wrappers
    PtrList<UOPstream> UOPstreamPtrs(Pstream::nProcs());

    // While there are particles to transfer
    while (true)
    {
        // Reset transfer buffers
        pBufs.clear();

        // Rewind existing streams
        forAll(UOPstreamPtrs, proci)
        {
            auto* osptr = UOPstreamPtrs.get(proci);
            if (osptr)
            {
                osptr->rewind();
            }
        }

        // Loop over all particles
        for (ParticleType& p : *this)
        {
            // Move the particle
            const bool keepParticle = p.move(cloud, td, trackTime);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                if (td.switchProcessor)
                {
                    #ifdef FULLDEBUG
                    if
                    (
                        !Pstream::parRun()
                     || !p.onBoundaryFace()
                     || procPatchNeighbours[p.patch()] < 0
                    )
                    {
                        FatalErrorInFunction
                            << "Switch processor flag is true when no parallel "
                            << "transfer is possible. This is a bug."
                            << exit(FatalError);
                    }
                    #endif

                    const label patchi = p.patch();

                    const label toProci =
                    (
                        refCast<const processorPolyPatch>(pbm[patchi])
                        .neighbProcNo()
                    );

                    // Get/create output stream
                    auto* osptr = UOPstreamPtrs.get(toProci);
                    if (!osptr)
                    {
                        osptr = new UOPstream(toProci, pBufs);
                        UOPstreamPtrs.set(toProci, osptr);
                    }

                    p.prepareForParallelTransfer();

                    // Tuple: (patchi particle)
                    (*osptr) << procPatchNeighbours[patchi] << p;

                    // Can now remove from my list
                    deleteParticle(p);
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        if (!Pstream::parRun())
        {
            break;
        }

        pBufs.finishedNeighbourSends(neighbourProcs);

        if (!returnReduceOr(pBufs.hasRecvData()))
        {
            // No parcels to transfer
            break;
        }

        // Retrieve from receive buffers
        for (const label proci : neighbourProcs)
        {
            if (pBufs.recvDataCount(proci))
            {
                UIPstream is(proci, pBufs);

                // Read out each (patchi particle) tuple
                while (!is.eof())
                {
                    label patchi = pTraits<label>(is);
                    auto* newp = new ParticleType(polyMesh_, is);

                    // The real patch index
                    patchi = procPatches[patchi];

                    (*newp).correctAfterParallelTransfer(patchi, td);
                    addParticle(newp);
                }
            }
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::autoMap(const mapPolyMesh& mapper)
{
    if (!globalPositionsPtr_)
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "Cloud::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    // Reset stored data that relies on the mesh
    //    polyMesh_.clearCellTree();
    cellWallFacesPtr_.clear();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    (void)polyMesh_.tetBasePtIs();
    (void)polyMesh_.oldCellCentres();

    const vectorField& positions = globalPositionsPtr_();

    label i = 0;
    for (ParticleType& p : *this)
    {
        p.autoMap(positions[i], mapper);
        ++i;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writePositions() const
{
    OFstream os
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    for (const ParticleType& p : *this)
    {
        const point position(p.position());
        os  << "v "
            << position.x() << ' '
            << position.y() << ' '
            << position.z() << nl;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::storeGlobalPositions() const
{
    // Store the global positions for later use by autoMap. It would be
    // preferable not to need this. If the mapPolyMesh object passed to autoMap
    // had a copy of the old mesh then the global positions could be recovered
    // within autoMap, and this pre-processing would not be necessary.

    globalPositionsPtr_.reset(new vectorField(this->size()));
    vectorField& positions = globalPositionsPtr_();

    label i = 0;
    for (const ParticleType& p : *this)
    {
        positions[i] = p.position();
        ++i;
    }
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
