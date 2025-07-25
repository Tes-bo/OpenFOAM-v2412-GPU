/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "enrichedPatch.H"
#include "OFstream.H"
#include "meshTools.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(enrichedPatch, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcMeshPoints() const
{
    if (meshPointsPtr_)
    {
        FatalErrorInFunction
            << "Mesh points already calculated."
            << abort(FatalError);
    }

    meshPointsPtr_.reset(new labelList(pointMap().sortedToc()));
}


void Foam::enrichedPatch::calcLocalFaces() const
{
    if (localFacesPtr_)
    {
        FatalErrorInFunction
            << "Local faces already calculated."
            << abort(FatalError);
    }

    // Invert mesh points and renumber faces using it
    const labelList& mp = meshPoints();

    Map<label> mpLookup(invertToMap(mp));


    // Create local faces.
    // Copy original faces and overwrite vertices after

    const faceList& faces = enrichedFaces();

    localFacesPtr_.reset(new faceList(faces));
    auto& locFaces = *localFacesPtr_;

    for (face& f : locFaces)
    {
        inplaceRenumber(mpLookup, f);
    }
}


void Foam::enrichedPatch::calcLocalPoints() const
{
    if (localPointsPtr_)
    {
        FatalErrorInFunction
            << "Local points already calculated."
            << abort(FatalError);
    }

    const labelList& mp = meshPoints();

    localPointsPtr_.reset(new pointField(mp.size()));
    auto& locPoints = *localPointsPtr_;

    forAll(locPoints, i)
    {
        locPoints[i] = *(pointMap().cfind(mp[i]));
    }
}


void Foam::enrichedPatch::clearOut()
{
    enrichedFacesPtr_.reset(nullptr);

    meshPointsPtr_.reset(nullptr);
    localFacesPtr_.reset(nullptr);
    localPointsPtr_.reset(nullptr);
    pointPointsPtr_.reset(nullptr);
    masterPointFacesPtr_.reset(nullptr);

    clearCutFaces();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enrichedPatch::enrichedPatch
(
    const primitiveFacePatch& masterPatch,
    const primitiveFacePatch& slavePatch,
    const labelUList& slavePointPointHits,
    const labelUList& slavePointEdgeHits,
    const UList<objectHit>& slavePointFaceHits
)
:
    masterPatch_(masterPatch),
    slavePatch_(slavePatch),
    pointMap_
    (
        masterPatch_.meshPoints().size()
      + slavePatch_.meshPoints().size()
    ),
    pointMapComplete_(false),
    pointMergeMap_(2*slavePatch_.meshPoints().size()),
    slavePointPointHits_(slavePointPointHits),
    slavePointEdgeHits_(slavePointEdgeHits),
    slavePointFaceHits_(slavePointFaceHits),
    enrichedFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    localFacesPtr_(nullptr),
    localPointsPtr_(nullptr),
    pointPointsPtr_(nullptr),
    masterPointFacesPtr_(nullptr),
    cutFacesPtr_(nullptr),
    cutFaceMasterPtr_(nullptr),
    cutFaceSlavePtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::enrichedPatch::meshPoints() const
{
    if (!meshPointsPtr_)
    {
        calcMeshPoints();
    }

    return *meshPointsPtr_;
}


const Foam::faceList& Foam::enrichedPatch::localFaces() const
{
    if (!localFacesPtr_)
    {
        calcLocalFaces();
    }

    return *localFacesPtr_;
}


const Foam::pointField& Foam::enrichedPatch::localPoints() const
{
    if (!localPointsPtr_)
    {
        calcLocalPoints();
    }

    return *localPointsPtr_;
}


const Foam::labelListList& Foam::enrichedPatch::pointPoints() const
{
    if (!pointPointsPtr_)
    {
        calcPointPoints();
    }

    return *pointPointsPtr_;
}


bool Foam::enrichedPatch::checkSupport() const
{
    const faceList& faces = enrichedFaces();

    bool error = false;

    forAll(faces, facei)
    {
        for (const label pointi : faces[facei])
        {
            if (!pointMap().found(pointi))
            {
                WarningInFunction
                    << "Point " << pointi << " of face " << facei
                    << " global point index: " << pointi
                    << " not supported in point map.  This is not allowed."
                    << endl;

                error = true;
            }
        }
    }

    return error;
}


void Foam::enrichedPatch::writeOBJ(const fileName& fName) const
{
    OFstream str(fName);

    meshTools::writeOBJ(str, localPoints());

    const faceList& faces = localFaces();

    for (const face& f : faces)
    {
        str << 'f';
        for (const label fp : f)
        {
            str << ' ' << fp+1;
        }
        str << nl;
    }
}


// ************************************************************************* //
