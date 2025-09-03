/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "fvFieldDecomposer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::fvFieldDecomposer::verbose_ = 1;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const labelUList& addressingSlice,
    const label addressingOffset
)
:
    directAddressing_(addressingSlice)
{
    forAll(directAddressing_, i)
    {
        // Subtract one to align addressing.
        directAddressing_[i] -= addressingOffset + 1;
    }
}


Foam::fvFieldDecomposer::processorVolPatchFieldDecomposer::
processorVolPatchFieldDecomposer
(
    const labelUList& owner,  // == mesh.faceOwner()
    const labelUList& neigh,  // == mesh.faceNeighbour()
    const labelUList& addressingSlice
)
:
    directAddressing_(addressingSlice.size())
{
    forAll(directAddressing_, i)
    {
        // Subtract one to align addressing.
        label ai = mag(addressingSlice[i]) - 1;

        if (ai < neigh.size())
        {
            // This is a regular face. it has been an internal face
            // of the original mesh and now it has become a face
            // on the parallel boundary.
            // Give face the value of the neighbour.

            if (addressingSlice[i] >= 0)
            {
                // I have the owner so use the neighbour value
                directAddressing_[i] = neigh[ai];
            }
            else
            {
                directAddressing_[i] = owner[ai];
            }
        }
        else
        {
            // This is a face that used to be on a cyclic boundary
            // but has now become a parallel patch face. I cannot
            // do the interpolation properly (I would need to look
            // up the different (face) list of data), so I will
            // just grab the value from the owner cell

            directAddressing_[i] = owner[ai];
        }
    }
}


Foam::fvFieldDecomposer::processorVolPatchFieldDecomposer::
processorVolPatchFieldDecomposer
(
    const polyMesh& mesh,
    const labelUList& addressingSlice
)
:
    processorVolPatchFieldDecomposer
    (
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        addressingSlice
    )
{}


Foam::fvFieldDecomposer::processorSurfacePatchFieldDecomposer::
processorSurfacePatchFieldDecomposer
(
    const labelUList& addressingSlice
)
:
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    forAll(addressing_, i)
    {
        addressing_[i].resize(1);
        weights_[i].resize(1);

        addressing_[i][0] = mag(addressingSlice[i]) - 1;
        weights_[i][0] = 1;
    }
}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const Foam::zero,
    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const labelList& boundaryAddressing
)
:
    procMesh_(procMesh),
    faceAddressing_(faceAddressing),
    cellAddressing_(cellAddressing),
    boundaryAddressing_(boundaryAddressing),
    // Mappers
    patchFieldDecomposerPtrs_(),
    processorVolPatchFieldDecomposerPtrs_(),
    processorSurfacePatchFieldDecomposerPtrs_(),
    faceSign_()
{}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const fvMesh& completeMesh,
    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const labelList& boundaryAddressing
)
:
    fvFieldDecomposer
    (
        zero{},
        procMesh,
        faceAddressing,
        cellAddressing,
        boundaryAddressing
    )
{
    reset(completeMesh);
}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const List<labelRange>& boundaryRanges,
    const labelUList& faceOwner,
    const labelUList& faceNeighbour,

    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const labelList& boundaryAddressing
)
:
    fvFieldDecomposer
    (
        zero{},
        procMesh,
        faceAddressing,
        cellAddressing,
        boundaryAddressing
    )
{
    reset(boundaryRanges, faceOwner, faceNeighbour);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvFieldDecomposer::empty() const
{
    return patchFieldDecomposerPtrs_.empty();
}


void Foam::fvFieldDecomposer::clear()
{
    patchFieldDecomposerPtrs_.clear();
    processorVolPatchFieldDecomposerPtrs_.clear();
    processorSurfacePatchFieldDecomposerPtrs_.clear();
    faceSign_.clear();
}


void Foam::fvFieldDecomposer::reset
(
    const List<labelRange>& boundaryRanges,
    const labelUList& faceOwner,
    const labelUList& faceNeighbour
)
{
    clear();
    const label nMappers = procMesh_.boundary().size();
    patchFieldDecomposerPtrs_.resize(nMappers);
    processorVolPatchFieldDecomposerPtrs_.resize(nMappers);
    processorSurfacePatchFieldDecomposerPtrs_.resize(nMappers);
    faceSign_.resize(nMappers);

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];
        const fvPatch& fvp = procMesh_.boundary()[patchi];
        const labelSubList localPatchSlice(fvp.patchSlice(faceAddressing_));

        if
        (
            oldPatchi >= 0
        && !isA<processorLduInterface>(procMesh_.boundary()[patchi])
        )
        {
            patchFieldDecomposerPtrs_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    localPatchSlice,
                    boundaryRanges[oldPatchi].start()
                )
            );
        }
        else
        {
            processorVolPatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorVolPatchFieldDecomposer
                (
                    faceOwner,
                    faceNeighbour,
                    localPatchSlice
                )
            );

            processorSurfacePatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorSurfacePatchFieldDecomposer
                (
                    static_cast<const labelUList&>(localPatchSlice)
                )
            );

            faceSign_.set
            (
                patchi,
                new scalarField(localPatchSlice.size())
            );

            {
                scalarField& s = faceSign_[patchi];
                forAll(s, i)
                {
                    s[i] = sign(localPatchSlice[i]);
                }
            }
        }
    }
}


void Foam::fvFieldDecomposer::reset(const fvMesh& completeMesh)
{
    clear();
    const label nMappers = procMesh_.boundary().size();
    patchFieldDecomposerPtrs_.resize(nMappers);
    processorVolPatchFieldDecomposerPtrs_.resize(nMappers);
    processorSurfacePatchFieldDecomposerPtrs_.resize(nMappers);
    faceSign_.resize(nMappers);

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];
        const fvPatch& fvp = procMesh_.boundary()[patchi];
        const labelSubList localPatchSlice(fvp.patchSlice(faceAddressing_));

        if
        (
            oldPatchi >= 0
        && !isA<processorLduInterface>(procMesh_.boundary()[patchi])
        )
        {
            patchFieldDecomposerPtrs_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    localPatchSlice,
                    completeMesh.boundaryMesh()[oldPatchi].start()
                )
            );
        }
        else
        {
            processorVolPatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorVolPatchFieldDecomposer
                (
                    completeMesh,
                    localPatchSlice
                )
            );

            processorSurfacePatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorSurfacePatchFieldDecomposer
                (
                    static_cast<const labelUList&>(localPatchSlice)
                )
            );

            faceSign_.set
            (
                patchi,
                new scalarField(localPatchSlice.size())
            );

            {
                scalarField& s = faceSign_[patchi];
                forAll(s, i)
                {
                    s[i] = sign(localPatchSlice[i]);
                }
            }
        }
    }
}


// ************************************************************************* //
