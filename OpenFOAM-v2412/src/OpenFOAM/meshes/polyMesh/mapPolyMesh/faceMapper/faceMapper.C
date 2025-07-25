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

#include "faceMapper.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceMapper::calcAddressing() const
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
            << "Addressing already calculated."
            << abort(FatalError);
    }

    if (direct())
    {
        // Direct addressing, no weights

        // Restrict addressing list to contain only live faces
        directAddrPtr_ = std::make_unique<labelList>
        (
            labelList::subList(mpm_.faceMap(), mapperLen_)
        );
        auto& directAddr = *directAddrPtr_;

        insertedObjectsPtr_ = std::make_unique<labelList>();
        auto& inserted = *insertedObjectsPtr_;

        // The nInsertedObjects_ already counted in the constructor
        if (nInsertedObjects_)
        {
            inserted.resize(nInsertedObjects_);

            label nInserted = 0;
            forAll(directAddr, i)
            {
                if (directAddr[i] < 0)
                {
                    // Found inserted
                    directAddr[i] = 0;
                    inserted[nInserted] = i;
                    ++nInserted;

                    // TBD: check (nInsertedObjects_ < nInserted)?
                    #ifdef FULLDEBUG
                    if (nInsertedObjects_ < nInserted)
                    {
                        FatalErrorInFunction
                            << "Unexpected insert of more than "
                            << nInsertedObjects_ << " items\n"
                            << abort(FatalError);
                    }
                    #endif
                }
            }
            // TBD: check (nInserted < nInsertedObjects_)?
            #ifdef FULLDEBUG
            if (nInserted < nInsertedObjects_)
            {
                WarningInFunction
                    << "Found " << nInserted << " instead of "
                    << nInsertedObjects_ << " items to insert\n";
            }
            #endif
            // The resize should be unnecessary
            inserted.resize(nInserted);
        }
    }
    else
    {
        // Interpolative addressing

        interpAddrPtr_ = std::make_unique<labelListList>(mapperLen_);
        auto& addr = *interpAddrPtr_;

        weightsPtr_ = std::make_unique<scalarListList>(mapperLen_);
        auto& wght = *weightsPtr_;


        // Set the addressing and uniform weight
        const auto setAddrWeights = [&]
        (
            const List<objectMap>& maps,
            const char * const nameOfMap
        )
        {
            for (const objectMap& map : maps)
            {
                // Get index, addressing
                const label facei = map.index();
                const labelList& mo = map.masterObjects();
                if (mo.empty()) continue;  // safety

                if (addr[facei].size())
                {
                    FatalErrorInFunction
                        << "Master face " << facei
                        << " already mapped, cannot apply "
                        << nameOfMap
                        << flatOutput(mo) << abort(FatalError);
                }

                // Map from masters, uniform weights
                addr[facei] = mo;
                wght[facei] = scalarList(mo.size(), 1.0/mo.size());
            }
        };


        setAddrWeights(mpm_.facesFromPointsMap(), "point faces");
        setAddrWeights(mpm_.facesFromEdgesMap(), "edge faces");
        setAddrWeights(mpm_.facesFromFacesMap(), "face faces");


        // Do mapped faces.
        // - may already have been set, so check if addressing still empty().

        {
            const labelList& map = mpm_.faceMap();

            // NB: faceMap can be longer than nFaces()
            for (label facei = 0; facei < mapperLen_; ++facei)
            {
                const label mappedi = map[facei];

                if (mappedi >= 0 && addr[facei].empty())
                {
                    // Mapped from a single face
                    addr[facei].resize(1, mappedi);
                    wght[facei].resize(1, 1.0);
                }
            }
        }


        // Grab inserted faces (for them the size of addressing is still zero)

        insertedObjectsPtr_ = std::make_unique<labelList>();
        auto& inserted = *insertedObjectsPtr_;

        // The nInsertedObjects_ already counted in the constructor
        if (nInsertedObjects_)
        {
            inserted.resize(nInsertedObjects_);

            label nInserted = 0;
            forAll(addr, i)
            {
                if (addr[i].empty())
                {
                    // Mapped from dummy face 0
                    addr[i].resize(1, 0);
                    wght[i].resize(1, 1.0);

                    inserted[nInserted] = i;
                    ++nInserted;

                    // TBD: check (nInsertedObjects_ < nInserted)?
                    #ifdef FULLDEBUG
                    if (nInsertedObjects_ < nInserted)
                    {
                        FatalErrorInFunction
                            << "Unexpected insert of more than "
                            << nInsertedObjects_ << " items\n"
                            << abort(FatalError);
                    }
                    #endif
                }
            }
            // TBD: check (nInserted < nInsertedObjects_)?
            #ifdef FULLDEBUG
            if (nInserted < nInsertedObjects_)
            {
                WarningInFunction
                    << "Found " << nInserted << " instead of "
                    << nInsertedObjects_ << " items to insert\n";
            }
            #endif
            // The resize should be unnecessary
            inserted.resize(nInserted);
        }
    }
}


// void Foam::faceMapper::clearOut()
// {
//     directAddrPtr_.reset(nullptr);
//     interpAddrPtr_.reset(nullptr);
//     weightsPtr_.reset(nullptr);
//     insertedObjectsPtr_.reset(nullptr);
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceMapper::faceMapper(const mapPolyMesh& mpm)
:
    mpm_(mpm),
    mapperLen_(mpm.mesh().nFaces()),
    nInsertedObjects_(0),
    direct_
    (
        // Mapping without interpolation?
        mpm.facesFromPointsMap().empty()
     && mpm.facesFromEdgesMap().empty()
     && mpm.facesFromFacesMap().empty()
    )
{
    const auto& directMap = mpm_.faceMap();

    if (!mapperLen_)
    {
        // Empty mesh
        direct_ = true;
        nInsertedObjects_ = 0;
    }
    else if (direct_)
    {
        // Number of inserted faces (-ve values)
        nInsertedObjects_ = std::count_if
        (
            directMap.cbegin(),
            directMap.cbegin(mapperLen_),
            [](label i) { return (i < 0); }
        );
    }
    else
    {
        // Check if there are inserted faces with no owner
        // (check all lists)

        bitSet unmapped(mapperLen_, true);

        unmapped.unset(directMap);  // direct mapped

        for (const auto& map : mpm_.facesFromPointsMap())
        {
            if (!map.empty()) unmapped.unset(map.index());
        }

        for (const auto& map : mpm_.facesFromEdgesMap())
        {
            if (!map.empty()) unmapped.unset(map.index());
        }

        for (const auto& map : mpm_.facesFromFacesMap())
        {
            if (!map.empty()) unmapped.unset(map.index());
        }

        nInsertedObjects_ = label(unmapped.count());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceMapper::~faceMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::faceMapper::size() const
{
    return mapperLen_;
}


Foam::label Foam::faceMapper::sizeBeforeMapping() const
{
    return mpm_.nOldFaces();
}


Foam::label Foam::faceMapper::internalSizeBeforeMapping() const
{
    return mpm_.nOldInternalFaces();
}


const Foam::labelUList& Foam::faceMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorInFunction
            << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!insertedObjects())
    {
        // No inserted faces.  Re-use faceMap
        return mpm_.faceMap();
    }
    else
    {
        if (!directAddrPtr_)
        {
            calcAddressing();
        }

        return *directAddrPtr_;
    }
}


const Foam::labelListList& Foam::faceMapper::addressing() const
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


const Foam::scalarListList& Foam::faceMapper::weights() const
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


const Foam::labelList& Foam::faceMapper::insertedObjectLabels() const
{
    if (!insertedObjectsPtr_)
    {
        if (!nInsertedObjects_)
        {
            // No inserted objects will be created
            return labelList::null();
        }

        calcAddressing();
    }

    return *insertedObjectsPtr_;
}


const Foam::labelHashSet& Foam::faceMapper::flipFaceFlux() const
{
    return mpm_.flipFaceFlux();
}


Foam::label Foam::faceMapper::nOldInternalFaces() const
{
    return mpm_.nOldInternalFaces();
}


const Foam::labelList& Foam::faceMapper::oldPatchStarts() const
{
    return mpm_.oldPatchStarts();
}


const Foam::labelList& Foam::faceMapper::oldPatchSizes() const
{
    return mpm_.oldPatchSizes();
}


// ************************************************************************* //
