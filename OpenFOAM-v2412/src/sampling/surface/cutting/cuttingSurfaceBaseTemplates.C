/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "volFields.H"
#include "edgeHashes.H"
#include "HashOps.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class EdgeOrientIntersect, class EdgeAlphaIntersect>
void Foam::cuttingSurfaceBase::walkCellCuts
(
    const primitiveMesh& mesh,
    const bitSet& cellCuts,
    const EdgeOrientIntersect& edgeOrientIntersect,
    const EdgeAlphaIntersect&  edgeAlphaIntersect,
    const bool triangulate,
    label nFaceCuts
)
{
    // Information required from the mesh
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const pointField& points = mesh.points();

    // Dynamic lists to handle triangulation and/or missed cuts etc
    const label nCellCuts = cellCuts.count();

    DynamicList<point> dynCutPoints(4*nCellCuts);
    DynamicList<face>  dynCutFaces(4*nCellCuts);
    DynamicList<label> dynCutCells(nCellCuts);

    // No nFaceCuts provided? Use a reasonable estimate
    if (!nFaceCuts)
    {
        nFaceCuts = 4*nCellCuts;
    }

    // Edge to pointId mapping
    EdgeMap<label> handledEdges(4*nFaceCuts);

    // Local scratch space for face vertices
    DynamicList<label> localFaceLoop(16);

    // Local scratch space for edge to pointId
    EdgeMap<label> localEdges;

    // Local scratch space for edge to faceId
    EdgeMap<edge>  localFaces;

    // Avoid duplicates for cuts exactly through a mesh point.
    // No other way to distinguish them, since there is no single edge
    // that "owns" the point.
    Map<label> endPoints;

    // Hash of faces (face points) that are exactly on a cell face
    HashSet<labelList> onCellFace;


    // Failure handling

    // Cells where walking failed (concave, degenerate, ...)
    labelHashSet failedCells;

    // To unwind insertion of end-point cutting
    labelHashSet localEndPoints;

    // Our recovery point on failure
    label unwindPoint = 0;

    // Cleanup routine for failed cell cut:
    const auto unwindWalk =
        [&](const label failedCellId = -1) -> void
        {
            // Discard points introduced
            dynCutPoints.resize(unwindPoint);

            // Discard end-point cuts
            endPoints.erase(localEndPoints);

            // Optionally record the failedCellId
            if (failedCellId != -1)
            {
                failedCells.insert(failedCellId);
            }
        };


    // Loop over cells that are cut
    for (const label celli : cellCuts)
    {
        const cell& cFace = cells[celli];

        // Reset local scratch
        localEdges.clear();
        localFaces.clear();
        localFaceLoop.clear();
        localEndPoints.clear();

        unwindPoint = dynCutPoints.size();


        // Classification for all the points cut - see intersectsFace() above
        // for more detail
        unsigned pointCutType = 0u;

        for (const label facei : cFace)
        {
            const face& f = faces[facei];

            forAll(f, fp)
            {
                edge e(f.edge(fp));

                // Action #1: Orient edge (+ve gradient) and detect intersect
                if (!edgeOrientIntersect(e))
                {
                    continue;
                }

                // Record the edge/face combination for the edge cut.
                // NB: the second operation is edge::insert() which places
                // facei in the unoccupied 'slot'
                localFaces(e).insert(facei);

                // Already handled cut point in this inner-loop?
                if (localEdges.found(e))
                {
                    // No new edge cut required
                    continue;
                }

                // Already handled cut point in the outer-loop?
                label cutPointId = handledEdges.lookup(e, -1);

                if (cutPointId >= 0)
                {
                    // Copy existing edge cut-point index
                    localEdges.insert(e, cutPointId);
                    continue;
                }

                // Expected id for the cut point
                cutPointId = dynCutPoints.size();

                const point& p0 = points[e.first()];
                const point& p1 = points[e.last()];

                // Action #2: edge cut alpha
                const scalar alpha = edgeAlphaIntersect(e);

                if (alpha < SMALL)
                {
                    pointCutType |= 0x1; // Cut at 0 (first)

                    const label endp = e.first();

                    if (endPoints.insert(endp, cutPointId))
                    {
                        localEndPoints.insert(endp);
                        dynCutPoints.append(p0);
                    }
                    else
                    {
                        cutPointId = endPoints[endp];
                    }
                }
                else if (alpha >= (1.0 - SMALL))
                {
                    pointCutType |= 0x2; // Cut at 1 (last)

                    const label endp = e.last();

                    if (endPoints.insert(endp, cutPointId))
                    {
                        localEndPoints.insert(endp);
                        dynCutPoints.append(p1);
                    }
                    else
                    {
                        cutPointId = endPoints[endp];
                    }
                }
                else
                {
                    pointCutType |= 0x4; // Cut between

                    dynCutPoints.append((1-alpha)*p0 + alpha*p1);
                }

                // Introduce new edge cut point
                localEdges.insert(e, cutPointId);
            }
        }


        // The keys of localEdges, localFaces are now identical.
        // The size of either should represent the number of points for
        // the resulting face loop.

        const label nTargetLoop = localFaces.size();

        if (nTargetLoop < 3)
        {
            unwindWalk(celli);
            continue;
        }


        // Handling cuts between two cells
        // After the previous edgeIntersectAndOrient call, the edge is oriented
        // according to the gradient.
        // If we only ever cut at the same edge end we know that we have
        // a cut coinciding with a cell face.

        if (pointCutType == 1 || pointCutType == 2)
        {
            // Hash the face-points to avoid duplicate faces
            if (!onCellFace.insert(HashTableOps::values(localEdges, true)))
            {
                DebugInfo
                    <<"skip duplicate on-place cut for cell " << celli
                    << " type (" << pointCutType << ")" << endl;

                // A duplicate is not failure
                unwindWalk();
                continue;
            }
        }


        // Start somewhere.

        // Since the local edges are oriented according to the gradient,
        // they can also be used to determine the correct face orientation.

        const edge refEdge = localFaces.begin().key();
        label nextFace = localFaces.begin().val()[0];

        DebugInfo
            << "search face " <<  nextFace << " IN " <<  localEdges << endl;

        for
        (
            label loopi = 0;
            localFaces.size() && loopi < 2*nTargetLoop;
            ++loopi
        )
        {
            bool ok = false;

            forAllIters(localFaces, iter)
            {
                DebugInfo
                    << "lookup " << nextFace << " in " << iter.val() << nl;

                // Find local index (0,1) or -1 on failure
                const label got = iter.val().which(nextFace);

                if (got != -1)
                {
                    ok = true;

                    // The other face
                    nextFace = iter.val()[(got?0:1)];

                    // The edge -> cut point
                    localFaceLoop.append(localEdges[iter.key()]);

                    DebugInfo
                        <<" faces " << iter.val()
                        << " point " << localFaceLoop.last()
                        << " edge=" << iter.key() << " nextFace=" << nextFace
                        << nl;

                    // Done this connection
                    localFaces.erase(iter);
                    break;
                }
            }

            if (!ok)
            {
                break;
            }
        }


        // Could also check if localFaces is now empty.

        if (nTargetLoop != localFaceLoop.size())
        {
            DebugInfo
                << "Warn expected " << nTargetLoop << " but got "
                << localFaceLoop.size() << endl;

            unwindWalk(celli);
            continue;
        }


        // Success
        handledEdges += localEdges;

        face f(localFaceLoop);

        // Orient face to point in the same direction as the edge gradient
        if ((f.areaNormal(dynCutPoints) & refEdge.vec(points)) < 0)
        {
            f.flip();
        }

        // The cut faces can be quite ugly, so optionally triangulate
        if (triangulate)
        {
            label nTri = f.triangles(dynCutPoints, dynCutFaces);
            while (nTri--)
            {
                dynCutCells.append(celli);
            }
        }
        else
        {
            dynCutFaces.append(f);
            dynCutCells.append(celli);
        }
    }

    if (failedCells.size())
    {
        WarningInFunction
            << "Failed cuts for " << failedCells.size() << " cells:" << nl
            << "    " << flatOutput(failedCells.sortedToc()) << nl
            << endl;
    }


    // No cuts? Then no need for any of this information
    if (dynCutCells.empty())
    {
        this->storedPoints().clear();
        this->storedFaces().clear();
        meshCells_.clear();
    }
    else
    {
        this->storedPoints().transfer(dynCutPoints);
        this->storedFaces().transfer(dynCutFaces);
        meshCells_.transfer(dynCutCells);
    }
}


// ************************************************************************* //
