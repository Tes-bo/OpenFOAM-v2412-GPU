/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "surfaceIntersection.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "edgeHashes.H"
#include "labelPairHashes.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// file-scope
// Write points in obj format
static void writeObjPoints(const UList<point>& pts, Ostream& os)
{
    for (const point& pt : pts)
    {
        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Get minimum length of all edges connected to point
Foam::scalar Foam::surfaceIntersection::minEdgeLen
(
    const triSurface& surf,
    const label pointi
)
{
    const labelList& pEdges = surf.pointEdges()[pointi];

    scalar minLen = GREAT;

    forAll(pEdges, pEdgeI)
    {
        const edge& e = surf.edges()[pEdges[pEdgeI]];

        minLen = min(minLen, e.mag(surf.localPoints()));
    }

    return minLen;
}


// Get edge between fp and fp+1 on facei.
Foam::label Foam::surfaceIntersection::getEdge
(
    const triSurface& surf,
    const label facei,
    const label fp
)
{
    const edge faceEdge = surf.localFaces()[facei].edge(fp);

    const labelList& eLabels = surf.faceEdges()[facei];

    forAll(eLabels, eI)
    {
        const label edgeI = eLabels[eI];

        if (surf.edges()[edgeI] == faceEdge)
        {
            return edgeI;
        }
    }

    FatalErrorInFunction
        << "Problem:: Cannot find edge with vertices " << faceEdge
        << " in face " << facei
        << abort(FatalError);

    return -1;
}


// Given a map remove all consecutive duplicate elements.
void Foam::surfaceIntersection::removeDuplicates
(
    const labelList& map,
    labelList& elems
)
{
    bool hasDuplicate = false;

    label prevVertI = -1;

    forAll(elems, elemI)
    {
        label newVertI = map[elems[elemI]];

        if (newVertI == prevVertI)
        {
            hasDuplicate = true;

            break;
        }
        prevVertI = newVertI;
    }

    if (hasDuplicate)
    {
        // Create copy
        labelList oldElems(elems);

        label elemI = 0;

        // Insert first
        elems[elemI++] = map[oldElems[0]];

        for (label vertI = 1; vertI < oldElems.size(); vertI++)
        {
            // Insert others only if they differ from one before
            label newVertI = map[oldElems[vertI]];

            if (newVertI != elems.last())
            {
                elems[elemI++] = newVertI;
            }
        }
        elems.setSize(elemI);
    }
}


// Remove all duplicate and degenerate elements. Return unique elements and
// map from old to new.
Foam::edgeList Foam::surfaceIntersection::filterEdges
(
    const edgeList& edges,
    labelList& map
)
{
    edgeHashSet uniqueEdges(10*edges.size());

    edgeList newEdges(edges.size());

    map.setSize(edges.size());
    map = -1;

    label newEdgeI = 0;

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if ((e.start() != e.end()) && uniqueEdges.insert(e))
        {
            // Edge is non-degenerate and not yet seen.

            map[edgeI] = newEdgeI;

            newEdges[newEdgeI++] = e;
        }
    }

    newEdges.setSize(newEdgeI);

    return newEdges;
}


// Remove all duplicate elements.
Foam::labelList Foam::surfaceIntersection::filterLabels
(
    const labelList& elems,
    labelList& map
)
{
    labelHashSet uniqueElems(10*elems.size());

    labelList newElems(elems.size());

    map.setSize(elems.size());
    map = -1;

    label newElemI = 0;

    forAll(elems, elemI)
    {
        label elem = elems[elemI];

        if (uniqueElems.insert(elem))
        {
            // First time elem is seen

            map[elemI] = newElemI;

            newElems[newElemI++] = elem;
        }
    }

    newElems.setSize(newElemI);

    return newElems;
}


void Foam::surfaceIntersection::writeIntersectedEdges
(
    const triSurface& surf,
    const labelListList& edgeCutVerts,
    Ostream& os
) const
{
    // Dump all points (surface followed by cutPoints)
    const pointField& pts = surf.localPoints();

    writeObjPoints(pts, os);
    writeObjPoints(cutPoints(), os);

    forAll(edgeCutVerts, edgeI)
    {
        const labelList& extraVerts = edgeCutVerts[edgeI];

        if (extraVerts.size())
        {
            const edge& e = surf.edges()[edgeI];

            // Start of original edge to first extra point
            os  << "l " << e.start()+1 << ' '
                << extraVerts[0] + surf.nPoints() + 1 << nl;

            for (label i = 1; i < extraVerts.size(); i++)
            {
                os  << "l " << extraVerts[i-1] + surf.nPoints() + 1  << ' '
                    << extraVerts[i] + surf.nPoints() + 1 << nl;
            }

            os  << "l " << extraVerts.last() + surf.nPoints() + 1
                << ' ' << e.end()+1 << nl;
        }
    }
}


// Return 0 (p close to start), 1(close to end) or -1.
Foam::label Foam::surfaceIntersection::classify
(
    const scalar startTol,
    const scalar endTol,
    const point& p,
    const edge& e,
    const UList<point>& points
)
{
    if (mag(p - points[e.start()]) < startTol)
    {
        return 0;
    }
    else if (mag(p - points[e.end()]) < endTol)
    {
        return 1;
    }

    return -1;
}


// ************************************************************************* //
