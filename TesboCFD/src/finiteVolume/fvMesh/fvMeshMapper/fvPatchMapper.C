/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "fvPatchMapper.H"
#include "fvPatch.H"
#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "mapPolyMesh.H"
#include "faceMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvPatchMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpAddrPtr_
     || weightsPtr_
    )
    {
        FatalErrorInFunction
            << "Addressing already calculated"
            << abort(FatalError);
    }

    // Mapping
    const label oldPatchStart =
        faceMap_.oldPatchStarts()[patch_.index()];

    const label oldPatchEnd =
        oldPatchStart + faceMap_.oldPatchSizes()[patch_.index()];

    hasUnmapped_ = false;

    // Assemble the maps: slice to patch
    if (direct())
    {
        // Direct mapping - slice to size
        directAddrPtr_ = std::make_unique<labelList>
        (
            patch_.patchSlice
            (
                static_cast<const labelList&>(faceMap_.directAddressing())
            )
        );
        auto& addr = *directAddrPtr_;

        // Adjust mapping to manage hits into other patches and into
        // internal
        forAll(addr, facei)
        {
            if
            (
                addr[facei] >= oldPatchStart
             && addr[facei] < oldPatchEnd
            )
            {
                addr[facei] -= oldPatchStart;
            }
            else
            {
                //addr[facei] = 0;
                addr[facei] = -1;
                hasUnmapped_ = true;
            }
        }

        if (fvMesh::debug)
        {
            if (min(addr) < 0)
            {
                WarningInFunction
                    << "Unmapped entry in patch mapping for patch "
                    << patch_.index() << " named " << patch_.name()
                    << endl;
            }
        }
    }
    else
    {
        // Interpolative mapping
        interpAddrPtr_ = std::make_unique<labelListList>
        (
            patch_.patchSlice(faceMap_.addressing())
        );
        auto& addr = *interpAddrPtr_;

        weightsPtr_ = std::make_unique<scalarListList>
        (
            patch_.patchSlice(faceMap_.weights())
        );
        auto& wght = *weightsPtr_;

        // Adjust mapping to manage hits into other patches and into
        // internal

        forAll(addr, facei)
        {
            auto& curAddr = addr[facei];
            auto& curWght = wght[facei];

            if
            (
                min(curAddr) >= oldPatchStart
             && max(curAddr) < oldPatchEnd
            )
            {
                // No adjustment of weights, just subtract patch start
                forAll(curAddr, i)
                {
                    curAddr[i] -= oldPatchStart;
                }
            }
            else
            {
                // Need to recalculate weights to exclude hits into internal

                label nActive = 0;
                scalar sumWeight = 0;

                forAll(curAddr, i)
                {
                    if
                    (
                        curAddr[i] >= oldPatchStart
                     && curAddr[i] < oldPatchEnd
                    )
                    {
                        curAddr[nActive] = curAddr[i] - oldPatchStart;
                        curWght[nActive] = curWght[i];

                        sumWeight += curWght[i];
                        ++nActive;
                    }
                }

                // Reset addressing and weights
                curAddr.resize(nActive);
                curWght.resize(nActive);

                // Re-scale the weights
                if (nActive)
                {
                    for (auto& w : curWght)
                    {
                        w /= sumWeight;
                    }
                }
                else
                {
                    hasUnmapped_ = true;
                }
            }
        }

        if (fvMesh::debug)
        {
            forAll(addr, i)
            {
                if (min(addr[i]) < 0)
                {
                    FatalErrorInFunction
                        << "Error in patch mapping for patch "
                        << patch_.index() << " named " << patch_.name()
                        << abort(FatalError);
                }
            }
        }
    }
}


// void Foam::fvPatchMapper::clearOut()
// {
//     directAddrPtr_.reset(nullptr);
//     interpAddrPtr_.reset(nullptr);
//     weightsPtr_.reset(nullptr);
//     hasUnmapped_ = false;
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvPatchMapper::fvPatchMapper
(
    const fvPatch& patch,
    const faceMapper& faceMap
)
:
    patch_(patch),
    faceMap_(faceMap),
    sizeBeforeMapping_(faceMap.oldPatchSizes()[patch_.index()]),
    hasUnmapped_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvPatchMapper::~fvPatchMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::fvPatchMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorInFunction
            << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::fvPatchMapper::addressing() const
{
    if (direct())
    {
        FatalErrorInFunction
            << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpAddrPtr_)
    {
        calcAddressing();
    }

    return *interpAddrPtr_;
}


const Foam::scalarListList& Foam::fvPatchMapper::weights() const
{
    if (direct())
    {
        FatalErrorInFunction
            << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


// ************************************************************************* //
