/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::enrichedPatch::enrichedFaceRatio_ = 3;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcEnrichedFaces
(
    const labelListList& pointsIntoMasterEdges,
    const labelListList& pointsIntoSlaveEdges,
    const pointField& projectedSlavePoints
)
{
    if (enrichedFacesPtr_)
    {
        FatalErrorInFunction
            << "Enriched faces already calculated."
            << abort(FatalError);
    }

    // Create a list of enriched faces
    // Algorithm:
    // 1) Grab the original face and start from point zero.
    // 2) If the point has been merged away, grab the merge label;
    //    otherwise, keep the original label.
    // 3) Go to the next edge. Collect all the points to be added along
    //    the edge; order them in the necessary direction and insert onto the
    //    face.
    // 4) Grab the next point and return on step 2.
    enrichedFacesPtr_.reset
    (
        new faceList(masterPatch_.size() + slavePatch_.size())
    );
    auto& enrichedFaces = *enrichedFacesPtr_;

    label nEnrichedFaces = 0;

    const pointField& masterLocalPoints = masterPatch_.localPoints();
    const faceList& masterLocalFaces = masterPatch_.localFaces();
    const labelListList& masterFaceEdges = masterPatch_.faceEdges();

    const faceList& slaveLocalFaces = slavePatch_.localFaces();
    const labelListList& slaveFaceEdges = slavePatch_.faceEdges();

    // For correct functioning of the enrichedPatch class, the slave
    // faces need to be inserted first.  See comments in
    // enrichedPatch.H

    // Get reference to the point merge map
    const Map<label>& pmm = pointMergeMap();

    // Add slave faces into the enriched faces list

    forAll(slavePatch_, facei)
    {
        const face& oldFace = slavePatch_[facei];
        const face& oldLocalFace = slaveLocalFaces[facei];
        const labelList& curEdges = slaveFaceEdges[facei];

        // Info<< "old slave face[" << facei << "] " << oldFace << endl;

        DynamicList<label> newFace(oldFace.size()*enrichedFaceRatio_);

        // Note: The number of points and edges in a face is always identical
        // so both can be done is the same loop
        forAll(oldFace, i)
        {
            // Add the point.
            // Using the mapped point id if possible

            const label mappedPointi = pmm.lookup(oldFace[i], oldFace[i]);

            newFace.append(mappedPointi);

            // Add the projected point into the patch support
            pointMap().insert
            (
                mappedPointi,   // Global label of point
                projectedSlavePoints[oldLocalFace[i]] // Projected position
            );

            // Grab the edge points

            const labelList& pointsOnEdge =
                pointsIntoSlaveEdges[curEdges[i]];

            // Info<< "slave pointsOnEdge for "
            //     << curEdges[i] << ": " << pointsOnEdge
            //     << endl;

            // If there are no points on the edge, skip everything
            // If there is only one point, no need for sorting
            if (pointsOnEdge.size())
            {
                // Sort edge points in order
                scalarField weightOnEdge(pointsOnEdge.size());

                const point& startPoint =
                    projectedSlavePoints[oldLocalFace[i]];

                const point& endPoint =
                    projectedSlavePoints[oldLocalFace.nextLabel(i)];

                vector e = (endPoint - startPoint);

                const scalar magSqrE = magSqr(e);

                if (magSqrE > SMALL)
                {
                    e /= magSqrE;
                }
                else
                {
                    FatalErrorInFunction
                        << "Zero length edge in slave patch for face " << i
                        << ".  This is not allowed."
                        << abort(FatalError);
                }

                pointField positionOnEdge(pointsOnEdge.size());

                forAll(pointsOnEdge, edgePointi)
                {
                    positionOnEdge[edgePointi] =
                        pointMap()[pointsOnEdge[edgePointi]];

                    weightOnEdge[edgePointi] =
                    (
                        e
                      &
                        (
                            positionOnEdge[edgePointi]
                          - startPoint
                        )
                    );
                }

                if (debug)
                {
                    // Check weights: all new points should be on the edge
                    if (min(weightOnEdge) < 0 || max(weightOnEdge) > 1)
                    {
                        FatalErrorInFunction
                            << " not on the edge for edge " << curEdges[i]
                            << " of face " << facei << " in slave patch." << nl
                            << "Min weight: " << min(weightOnEdge)
                            << " Max weight: " << max(weightOnEdge)
                            << abort(FatalError);
                    }
                }

                // Go through the points and collect them based on
                // weights from lower to higher.  This gives the
                // correct order of points along the edge.
                forAll(weightOnEdge, passI)
                {
                    // Max weight can only be one, so the sorting is
                    // done by elimination.
                    label nextPoint = -1;
                    scalar dist = 2;

                    forAll(weightOnEdge, wI)
                    {
                        if (weightOnEdge[wI] < dist)
                        {
                            dist = weightOnEdge[wI];
                            nextPoint = wI;
                        }
                    }

                    // Insert the next point and reset its weight to exclude it
                    // from future picks
                    newFace.append(pointsOnEdge[nextPoint]);
                    weightOnEdge[nextPoint] = GREAT;

                    // Add the point into patch support
                    pointMap().insert
                    (
                        pointsOnEdge[nextPoint],
                        positionOnEdge[nextPoint]
                    );
                }
            }
        }

        // Info<< "New slave face[" << facei << "] "
        //     << flatOutput(newFace) << " was " << flatOutput(oldFace)
        //     << endl;

        // Add the new face to the list
        enrichedFaces[nEnrichedFaces].transfer(newFace);
        nEnrichedFaces++;
    }


    // Add master faces into the enriched faces list

    forAll(masterPatch_, facei)
    {
        const face& oldFace = masterPatch_[facei];
        const face& oldLocalFace = masterLocalFaces[facei];
        const labelList& curEdges = masterFaceEdges[facei];

        // Info<< "old master face[" << facei << "] " << oldFace << endl;

        DynamicList<label> newFace(oldFace.size()*enrichedFaceRatio_);

        // Note: The number of points and edges in a face is always identical
        // so both can be done is the same loop
        forAll(oldFace, i)
        {
            // Add the point.
            // Using the mapped point id if possible

            const label mappedPointi = pmm.lookup(oldFace[i], oldFace[i]);

            newFace.append(mappedPointi);

            // Add the point into patch support
            pointMap().insert
            (
                mappedPointi,   // Global label of point
                masterLocalPoints[oldLocalFace[i]]
            );

            // Grab the edge points

            const labelList& pointsOnEdge =
                pointsIntoMasterEdges[curEdges[i]];

            // If there are no points on the edge, skip everything
            // If there is only one point, no need for sorting
            if (pointsOnEdge.size())
            {
                // Sort edge points in order
                scalarField weightOnEdge(pointsOnEdge.size());

                const point& startPoint =
                    masterLocalPoints[oldLocalFace[i]];

                const point& endPoint =
                    masterLocalPoints[oldLocalFace.nextLabel(i)];

                vector e = (endPoint - startPoint);

                const scalar magSqrE = magSqr(e);

                if (magSqrE > SMALL)
                {
                    e /= magSqrE;
                }
                else
                {
                    FatalErrorInFunction
                        << "Zero length edge in master patch for face " << i
                        << ".  This is not allowed."
                        << abort(FatalError);
                }

                pointField positionOnEdge(pointsOnEdge.size());

                forAll(pointsOnEdge, edgePointi)
                {
                    positionOnEdge[edgePointi] =
                        pointMap()[pointsOnEdge[edgePointi]];

                    weightOnEdge[edgePointi] =
                    (
                        e
                      &
                        (
                            positionOnEdge[edgePointi] - startPoint
                        )
                    );
                }

                if (debug)
                {
                    // Check weights: all new points should be on the edge
                    if (min(weightOnEdge) < 0 || max(weightOnEdge) > 1)
                    {
                        FatalErrorInFunction
                            << " not on the edge for edge " << curEdges[i]
                            << " of face " << facei << " in master patch." << nl
                            << "Min weight: " << min(weightOnEdge)
                            << " Max weight: " << max(weightOnEdge)
                            << abort(FatalError);
                    }
                }

                // Go through the points and collect them based on
                // weights from lower to higher.  This gives the
                // correct order of points along the edge.
                forAll(weightOnEdge, passI)
                {
                    // Max weight can only be one, so the sorting is
                    // done by elimination.
                    label nextPoint = -1;
                    scalar dist = 2;

                    forAll(weightOnEdge, wI)
                    {
                        if (weightOnEdge[wI] < dist)
                        {
                            dist = weightOnEdge[wI];
                            nextPoint = wI;
                        }
                    }

                    // Insert the next point and reset its weight to exclude it
                    // from future picks
                    newFace.append(pointsOnEdge[nextPoint]);
                    weightOnEdge[nextPoint] = GREAT;

                    // Add the point into patch support
                    pointMap().insert
                    (
                        pointsOnEdge[nextPoint],
                        positionOnEdge[nextPoint]
                    );
                }
            }
        }

        // Info<< "New master face[" << facei << "] "
        //     << flatOutput(newFace) << " was " << flatOutput(oldFace)
        //     << endl;

        // Add the new face to the list
        enrichedFaces[nEnrichedFaces].transfer(newFace);
        nEnrichedFaces++;
    }

    // Check the support for the enriched patch
    if (debug)
    {
        if (!checkSupport())
        {
            Info<< "Enriched patch support OK. Slave faces: "
                << slavePatch_.size() << " Master faces: "
                << masterPatch_.size() << endl;
        }
        else
        {
            FatalErrorInFunction
                << "Error in enriched patch support"
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faceList& Foam::enrichedPatch::enrichedFaces() const
{
    if (!enrichedFacesPtr_)
    {
        FatalErrorInFunction
            << "void enrichedPatch::calcEnrichedFaces\n"
            << "(\n"
            << "    const labelListList& pointsIntoMasterEdges,\n"
            << "    const labelListList& pointsIntoSlaveEdges,\n"
            << "    const pointField& projectedSlavePoints\n"
            << ")"
            << " before trying to access faces."
            << abort(FatalError);
    }

    return *enrichedFacesPtr_;
}


// ************************************************************************* //
