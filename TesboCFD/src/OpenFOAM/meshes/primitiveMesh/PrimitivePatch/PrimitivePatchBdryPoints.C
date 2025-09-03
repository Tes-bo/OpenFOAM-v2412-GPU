/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "PrimitivePatch.H"
#include "HashSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcBdryPoints() const
{
    if (boundaryPointsPtr_)
    {
        // Error to recalculate if already allocated
        FatalErrorInFunction
            << "boundaryPoints already calculated"
            << abort(FatalError);
    }

    labelHashSet bp(0);

    if (hasEdges())
    {
        DebugInFunction
            << "Calculating boundary points from existing addressing"
            << nl;

        bp.reserve(2*nBoundaryEdges());

        for (const edge& e : boundaryEdges())
        {
            bp.insert(e.first());
            bp.insert(e.second());
        }
    }
    else
    {
        DebugInFunction
            << "Calculating boundary points with manual edge addressing"
            << nl;


        // Calculate manually.
        // Needs localFaces, but uses local hashes of the edges here
        // instead of forcing a full faceFaces/edgeFaces/faceEdges calculation

        // Get reference to localFaces
        const List<face_type>& locFcs = localFaces();

        // Guess the max number of edges/neighbours for a face
        label edgeCount = 0;
        for (const auto& f : locFcs)
        {
            edgeCount += f.nEdges();
        }

        // ie, EdgeMap<label> to keep counts
        HashTable<label, edge, Hash<edge>> knownEdges(2*edgeCount);

        for (const auto& f : locFcs)
        {
            const label numEdges = f.nEdges();

            for (label edgei = 0; edgei < numEdges; ++edgei)
            {
                ++ knownEdges(f.edge(edgei));
            }
        }

        edgeCount = 0;

        forAllConstIters(knownEdges, iter)
        {
            if (1 == iter.val())  // Singly connected edge
            {
                ++edgeCount;
            }
        }

        bp.reserve(2*edgeCount);

        forAllConstIters(knownEdges, iter)
        {
            const edge& e = iter.key();

            if (1 == iter.val())  // Singly connected edge
            {
                bp.insert(e.first());
                bp.insert(e.second());
            }
        }
    }

    boundaryPointsPtr_.reset(new labelList(bp.sortedToc()));
    DebugInfo << "    Finished." << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::Pair<typename Foam::PrimitivePatch<FaceList, PointField>::point_type>
Foam::PrimitivePatch<FaceList, PointField>::box() const
{
    Pair<point_type> bb
    (
        point_type::uniform(ROOTVGREAT),
        point_type::uniform(-ROOTVGREAT)
    );

    if (hasMeshPoints())
    {
        // Less looping if meshPoints() are already available
        for (const label pointi : meshPoints())
        {
            bb.first()  = min(bb.first(),  points_[pointi]);
            bb.second() = max(bb.second(), points_[pointi]);
        }
    }
    else
    {
        // Walk the points on each face
        for (const face_type& f : *this)
        {
            for (const label pointi : f)
            {
                bb.first()  = min(bb.first(),  points_[pointi]);
                bb.second() = max(bb.second(), points_[pointi]);
            }
        }
    }

    return bb;
}


template<class FaceList, class PointField>
Foam::scalar
Foam::PrimitivePatch<FaceList, PointField>::sphere(const label facei) const
{
    scalar radiusSqr = 0;

    const point_type& fc = this->faceCentres()[facei];

    for (const label fp : this->operator[](facei))
    {
        const scalar sqrDist = magSqr(fc - points_[fp]);
        if (radiusSqr < sqrDist)
        {
            radiusSqr = sqrDist;
        }
    }

    return radiusSqr;
}


// ************************************************************************* //
