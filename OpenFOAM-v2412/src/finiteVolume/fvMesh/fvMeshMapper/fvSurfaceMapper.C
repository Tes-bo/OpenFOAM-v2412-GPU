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

#include "fvSurfaceMapper.H"
#include "fvMesh.H"
#include "mapPolyMesh.H"
#include "faceMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvSurfaceMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpAddrPtr_
     || weightsPtr_
     || insertedObjectsPtr_
    )
    {
        FatalErrorInFunction
            << "Addressing already calculated"
            << abort(FatalError);
    }

    // Mapping

    const label oldNInternal = faceMap_.nOldInternalFaces();

    // Assemble the maps
    if (direct())
    {
        // Direct mapping - slice to size
        directAddrPtr_ = std::make_unique<labelList>
        (
            labelList::subList(faceMap_.directAddressing(), size())
        );
        auto& addr = *directAddrPtr_;

        // Adjust for creation of an internal face from a boundary face
        forAll(addr, facei)
        {
            if (addr[facei] > oldNInternal)
            {
                addr[facei] = 0;
            }
        }
    }
    else
    {
        // Interpolative mapping - slice to size
        interpAddrPtr_ = std::make_unique<labelListList>
        (
            labelListList::subList(faceMap_.addressing(), size())
        );
        auto& addr = *interpAddrPtr_;

        weightsPtr_ = std::make_unique<scalarListList>
        (
            scalarListList::subList(faceMap_.weights(), size())
        );
        auto& wght = *weightsPtr_;

        // Adjust for creation of an internal face from a boundary face
        forAll(addr, facei)
        {
            if (max(addr[facei]) >= oldNInternal)
            {
                addr[facei] = labelList(1, Foam::zero{});
                wght[facei] = scalarList(1, scalar(1));
            }
        }
    }

    // Inserted objects

    insertedObjectsPtr_ = std::make_unique<labelList>();
    auto& inserted = *insertedObjectsPtr_;

    // If there are, assemble the labels
    if (faceMap_.insertedObjects())
    {
        const labelList& insFaces = faceMap_.insertedObjectLabels();

        inserted.resize(insFaces.size());

        label count = 0;
        for (const label facei : insFaces)
        {
            // If the face is internal, keep it here
            if (facei < size())
            {
                inserted[count] = facei;
                ++count;
            }
        }

        inserted.resize(count);
    }
}


// void Foam::fvSurfaceMapper::clearOut()
// {
//     directAddrPtr_.reset(nullptr);
//     interpAddrPtr_.reset(nullptr);
//     weightsPtr_.reset(nullptr);
//     insertedObjectsPtr_.reset(nullptr);
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvSurfaceMapper::fvSurfaceMapper
(
    const fvMesh& mesh,
    const faceMapper& mapper
)
:
    mesh_(mesh),
    faceMap_(mapper)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvSurfaceMapper::~fvSurfaceMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::fvSurfaceMapper::directAddressing() const
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


const Foam::labelListList& Foam::fvSurfaceMapper::addressing() const
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


const Foam::scalarListList& Foam::fvSurfaceMapper::weights() const
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


const Foam::labelList& Foam::fvSurfaceMapper::insertedObjectLabels() const
{
    if (!insertedObjectsPtr_)
    {
        calcAddressing();
    }

    return *insertedObjectsPtr_;
}


// ************************************************************************* //
