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

#include "pointMapper.H"
#include "pointMesh.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMapper::calcAddressing() const
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

        directAddrPtr_ = std::make_unique<labelList>
        (
            // No retired points, so pointMap().size() == mapperLen_ anyhow
            labelList::subList(mpm_.pointMap(), mapperLen_)
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
                const label pointi = map.index();
                const labelList& mo = map.masterObjects();
                if (mo.empty()) continue;  // safety

                if (addr[pointi].size())
                {
                    FatalErrorInFunction
                        << "Master point " << pointi
                        << " already mapped, cannot apply "
                        << nameOfMap
                        << flatOutput(mo) << abort(FatalError);
                }

                // Map from masters, uniform weights
                addr[pointi] = mo;
                wght[pointi] = scalarList(mo.size(), 1.0/mo.size());
            }
        };


        // Points created from other points (i.e. points merged into it).

        setAddrWeights(mpm_.pointsFromPointsMap(), "point points");


        // Do mapped points.
        // - may already have been set, so check if addressing still empty().

        {
            const labelList& map = mpm_.pointMap();

            for (label pointi = 0; pointi < mapperLen_; ++pointi)
            {
                const label mappedi = map[pointi];

                if (mappedi >= 0 && addr[pointi].empty())
                {
                    // Mapped from a single point
                    addr[pointi].resize(1, mappedi);
                    wght[pointi].resize(1, 1.0);
                }
            }
        }

        // Grab inserted points (for them the size of addressing is still zero)

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
                    // Mapped from dummy point 0
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


// void Foam::pointMapper::clearOut()
// {
//     directAddrPtr_.reset(nullptr);
//     interpAddrPtr_.reset(nullptr);
//     weightsPtr_.reset(nullptr);
//     insertedObjectsPtr_.reset(nullptr);
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMapper::pointMapper(const pointMesh& pMesh, const mapPolyMesh& mpm)
:
    mpm_(mpm),
    mapperLen_(pMesh.size()),
    nInsertedObjects_(0),
    direct_
    (
        // Mapping without interpolation?
        mpm.pointsFromPointsMap().empty()
    )
{
    const auto& directMap = mpm_.pointMap();

    if (!mapperLen_)
    {
        // Empty mesh
        direct_ = true;
        nInsertedObjects_ = 0;
    }
    else if (direct_)
    {
        // Number of inserted points (-ve values)
        nInsertedObjects_ = std::count_if
        (
            directMap.cbegin(),
            directMap.cbegin(mapperLen_),
            [](label i) { return (i < 0); }
        );
    }
    else
    {
        // Check if there are inserted points with no owner
        // (check all lists)

        bitSet unmapped(mapperLen_, true);

        unmapped.unset(directMap);  // direct mapped

        for (const objectMap& map : mpm_.pointsFromPointsMap())
        {
            if (!map.empty()) unmapped.unset(map.index());
        }

        nInsertedObjects_ = label(unmapped.count());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMapper::~pointMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pointMapper::size() const
{
    // OR:  return mapperLen_;
    return mpm_.pointMap().size();
}


Foam::label Foam::pointMapper::sizeBeforeMapping() const
{
    return mpm_.nOldPoints();
}


const Foam::labelUList& Foam::pointMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorInFunction
            << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!insertedObjects())
    {
        // No inserted points.  Re-use pointMap
        return mpm_.pointMap();
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


const Foam::labelListList& Foam::pointMapper::addressing() const
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


const Foam::scalarListList& Foam::pointMapper::weights() const
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


const Foam::labelList& Foam::pointMapper::insertedObjectLabels() const
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


// ************************************************************************* //
