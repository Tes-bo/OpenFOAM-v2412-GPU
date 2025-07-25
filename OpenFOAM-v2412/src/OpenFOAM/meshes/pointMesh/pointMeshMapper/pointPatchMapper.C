/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "pointPatchMapper.H"
#include "pointPatch.H"
#include "mapPolyMesh.H"
#include "faceMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointPatchMapper::calcAddressing() const
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

    hasUnmapped_ = false;

    if (direct())
    {
        // Direct mapping.
        directAddrPtr_ = std::make_unique<labelList>
        (
            mpm_.patchPointMap()[patch_.index()]
        );
        auto& addr = *directAddrPtr_;

        forAll(addr, i)
        {
            if (addr[i] < 0)
            {
                hasUnmapped_ = true;
                break;
            }
        }
    }
    else
    {
        // Interpolative mapping.

        // NOTE: Incorrect:
        // FOR NOW only takes first patch point instead of averaging all
        // patch points. Problem is we don't know what points were in the patch
        // for points that were merged.

        interpAddrPtr_ = std::make_unique<labelListList>(size());
        auto& addr = *interpAddrPtr_;

        weightsPtr_ = std::make_unique<scalarListList>(addr.size());
        auto& wght = *weightsPtr_;

        const labelList& ppm = mpm_.patchPointMap()[patch_.index()];

        forAll(ppm, i)
        {
            if (ppm[i] >= 0)
            {
                addr[i].resize(1, ppm[i]);
                wght[i].resize(1, 1.0);
            }
            else
            {
                // Inserted point.
                ///// Map from point0 (arbitrary choice)
                //addr[i].resize(1, 0);
                //wght[i].resize(1, 1.0);
                hasUnmapped_ = true;
            }
        }
    }
}


// void Foam::pointPatchMapper::clearOut()
// {
//     directAddrPtr_.reset(nullptr);
//     interpAddrPtr_.reset(nullptr);
//     weightsPtr_.reset(nullptr);
//     hasUnmapped_ = false;
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointPatchMapper::pointPatchMapper
(
    const pointPatch& patch,
    const pointMapper& pointMap,
    const mapPolyMesh& mpm
)
:
    pointPatchFieldMapper(),
    patch_(patch),
    pointMapper_(pointMap),
    mpm_(mpm),
    sizeBeforeMapping_
    (
        patch_.index() < mpm_.oldPatchNMeshPoints().size()
      ? mpm_.oldPatchNMeshPoints()[patch_.index()]
      : 0
    ),
    hasUnmapped_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointPatchMapper::~pointPatchMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::pointPatchMapper::directAddressing() const
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


const Foam::labelListList& Foam::pointPatchMapper::addressing() const
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


const Foam::scalarListList& Foam::pointPatchMapper::weights() const
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
