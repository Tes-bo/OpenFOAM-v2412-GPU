/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "faAreaMapper.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faAreaMapper::calcAddressing() const
{
    if
    (
        newFaceLabelsPtr_
     || newFaceLabelsMapPtr_
     || directAddrPtr_
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

    const label oldNInternal = mpm_.nOldInternalFaces();

    hasUnmapped_ = false;

    // Calculate new face labels

    // Copy old face labels
    const labelList& oldFaces = mesh_.faceLabels();

    // Prepare a list of new face labels and (preliminary) addressing
    // Note: dimensioned to number of boundary faces of polyMesh
    newFaceLabelsPtr_ = std::make_unique<labelList>
    (
        mesh_.mesh().nBoundaryFaces(),
        -1
    );
    auto& newFaceLabels = *newFaceLabelsPtr_;

    newFaceLabelsMapPtr_ = std::make_unique<labelList>
    (
        mesh_.mesh().nBoundaryFaces(),
        -1
    );
    auto& newFaceLabelsMap = *newFaceLabelsMapPtr_;
    label nNewFaces = 0;

    Info<< "Old face list size: " << oldFaces.size()
        << " estimated new size " << newFaceLabels.size() << endl;

    // Get reverse face map
    const labelList& reverseFaceMap = mpm_.reverseFaceMap();

    // Pick up live old faces
    forAll(oldFaces, faceI)
    {
        if (reverseFaceMap[oldFaces[faceI]] > -1)
        {
            // Face is live, add it and record addressing
            newFaceLabels[nNewFaces] = reverseFaceMap[oldFaces[faceI]];
            newFaceLabelsMap[nNewFaces] = faceI;

            ++nNewFaces;
        }
    }

    // Assemble the maps
    if (direct())
    {
        Info<< "Direct"<< endl;
        // Direct mapping: no further faces to add.  Resize list
        newFaceLabels.resize(nNewFaces);

        directAddrPtr_ = std::make_unique<labelList>
        (
            labelList::subList(newFaceLabelsMap, nNewFaces)
        );
        auto& addr = *directAddrPtr_;

        // Adjust for creation of a boundary face from an internal face
        forAll(addr, facei)
        {
            if (addr[facei] < oldNInternal)
            {
                addr[facei] = 0;
            }
        }
    }
    else
    {
        // There are further faces to add.  Prepare interpolation addressing
        // and weights to full size
        interpAddrPtr_ = std::make_unique<labelListList>
        (
            newFaceLabels.size()
        );
        auto& addr = *interpAddrPtr_;

        weightsPtr_ = std::make_unique<scalarListList>(addr.size());
        auto& wght = *weightsPtr_;

        // Insert single addressing and weights
        for (label addrI = 0; addrI < nNewFaces; ++addrI)
        {
            addr[addrI].resize(1, newFaceLabelsMap[addrI]);
            wght[addrI].resize(1, 1.0);
        }

        // Pick up faces from points, edges and faces where the origin
        // Only map from faces which were previously in the faMesh, using
        // fast lookup

        // Set of faces previously in the mesh
        const labelHashSet oldFaceLookup(oldFaces);

        // Check if master objects are in faMesh
        DynamicList<label> validMo(128);

        const auto addCheckedObjects = [&](const List<objectMap>& maps)
        {
            for (const objectMap& map : maps)
            {
                // Get target index, addressing
                const label facei = map.index();
                const labelList& mo = map.masterObjects();
                if (mo.empty()) continue;  // safety

                validMo.clear();
                validMo.reserve(mo.size());

                for (const label obji : mo)
                {
                    if (oldFaceLookup.contains(obji))
                    {
                        validMo.push_back(obji);
                    }
                }

                if (validMo.size())
                {
                    // Some objects found: add face and interpolation to list
                    newFaceLabels[nNewFaces] = facei;

                    // No old face available
                    newFaceLabelsMap[nNewFaces] = -1;

                    // Map from masters, uniform weights
                    addr[nNewFaces] = validMo;
                    wght[nNewFaces] =
                        scalarList(validMo.size(), 1.0/validMo.size());

                    ++nNewFaces;
                }
            }
        };


        // Go through faces-from lists and add the ones where all
        // old face labels belonged to the faMesh

        {
            addCheckedObjects(mpm_.facesFromPointsMap());
            addCheckedObjects(mpm_.facesFromEdgesMap());
            addCheckedObjects(mpm_.facesFromFacesMap());
        }

        // All faces collected.  Reset sizes of lists
        newFaceLabels.resize(nNewFaces);
        newFaceLabelsMap.resize(nNewFaces);
        addr.resize(nNewFaces);
        wght.resize(nNewFaces);

        Info<< "addr: " << addr << nl
            << "wght: " << wght << endl;
    }

    // Inserted objects cannot appear in the new faMesh as they have no master
    // HJ, 10/Aug/2011
    insertedObjectsPtr_ = std::make_unique<labelList>();
}


// void Foam::faAreaMapper::clearOut()
// {
//     newFaceLabelsPtr_.reset(nullptr);
//     newFaceLabelsMapPtr_.reset(nullptr);
//
//     directAddrPtr_.reset(nullptr);
//     interpAddrPtr_.reset(nullptr);
//     weightsPtr_.reset(nullptr);
//
//     insertedObjectsPtr_.reset(nullptr);
//     hasUnmapped_ = false;
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faAreaMapper::faAreaMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    mpm_(mpm),
    sizeBeforeMapping_(mesh.nFaces()),
    direct_
    (
        // Mapping without interpolation?
        mpm.facesFromPointsMap().empty()
     && mpm.facesFromEdgesMap().empty()
     && mpm.facesFromFacesMap().empty()
    ),
    hasUnmapped_(false)
{
    // Inserted objects not supported: no master
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faAreaMapper::~faAreaMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::faAreaMapper::newFaceLabels() const
{
    if (!newFaceLabelsPtr_)
    {
        calcAddressing();
    }

    return *newFaceLabelsPtr_;
}


const Foam::labelList& Foam::faAreaMapper::newFaceLabelsMap() const
{
    if (!newFaceLabelsMapPtr_)
    {
        calcAddressing();
    }

    return *newFaceLabelsMapPtr_;
}


const Foam::labelUList& Foam::faAreaMapper::directAddressing() const
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


const Foam::labelListList& Foam::faAreaMapper::addressing() const
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


const Foam::scalarListList& Foam::faAreaMapper::weights() const
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


const Foam::labelList& Foam::faAreaMapper::insertedObjectLabels() const
{
    if (!insertedObjectsPtr_)
    {
        calcAddressing();
    }

    return *insertedObjectsPtr_;
}


// ************************************************************************* //
