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

Description
    For every point on the patch find the closest face on the target side.
    Return a target face label for each patch point

\*---------------------------------------------------------------------------*/

#include "boolList.H"
#include "pointHit.H"
#include "objectHit.H"
#include "bandCompression.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
template<class ToPatch>
Foam::List<Foam::objectHit>
Foam::PrimitivePatch<FaceList, PointField>::projectPoints
(
    const ToPatch& targetPatch,
    const Field
    <
        typename Foam::PrimitivePatch<FaceList, PointField>::point_type
    >& projectionDirection,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // The current patch is slave, i.e. it is being projected onto the target

    if (projectionDirection.size() != nPoints())
    {
        FatalErrorInFunction
            << "Projection direction field does not correspond to "
            << "patch points." << endl
            << "Size: " << projectionDirection.size()
            << " Number of points: " << nPoints()
            << abort(FatalError);
    }

    const labelList& slavePointOrder = localPointOrder();

    const labelList& slaveMeshPoints = meshPoints();

    // Result
    List<objectHit> result(nPoints());

    const labelListList& masterFaceFaces = targetPatch.faceFaces();

    const ToPatch& masterFaces = targetPatch;

    const Field<point_type>& masterPoints = targetPatch.points();

    // Estimate face centre of target side
    Field<point_type> masterFaceCentres(targetPatch.size());

    forAll(masterFaceCentres, facei)
    {
        masterFaceCentres[facei] =
            average(masterFaces[facei].points(masterPoints));
    }

    // Algorithm:
    // Loop through all points of the slave side. For every point find the
    // radius for the current contact face. If the contact point falls inside
    // the face and the radius is smaller than for all neighbouring faces,
    // the contact is found. If not, visit the neighbour closest to the
    // calculated contact point. If a single master face is visited more than
    // twice, initiate n-squared search.

    label curFace = 0;
    label nNSquaredSearches = 0;

    forAll(slavePointOrder, pointi)
    {
        // Pick up slave point and direction
        const label curLocalPointLabel = slavePointOrder[pointi];

        const point_type& curPoint =
            points_[slaveMeshPoints[curLocalPointLabel]];

        const point_type& curProjectionDir =
            projectionDirection[curLocalPointLabel];

        bool closer;

        boolList visitedTargetFace(targetPatch.size(), false);
        bool doNSquaredSearch = false;

        bool foundEligible = false;

        scalar sqrDistance = GREAT;

        // Force the full search for the first point to ensure good
        // starting face
        if (pointi == 0)
        {
            doNSquaredSearch = true;
        }
        else
        {
            do
            {
                closer = false;
                doNSquaredSearch = false;

                // Calculate intersection with curFace
                PointHit<point_type> curHit =
                    masterFaces[curFace].ray
                    (
                        curPoint,
                        curProjectionDir,
                        masterPoints,
                        alg,
                        dir
                    );

                visitedTargetFace[curFace] = true;

                if (curHit.hit())
                {
                    result[curLocalPointLabel] = objectHit(true, curFace);

                    break;
                }
                else
                {
                    // If a new miss is eligible, it is closer than
                    // any previous eligible miss (due to surface walk)

                    // Only grab the miss if it is eligible
                    if (curHit.eligibleMiss())
                    {
                        foundEligible = true;
                        result[curLocalPointLabel] = objectHit(false, curFace);
                    }

                    // Find the next likely face for intersection

                    // Calculate the miss point on the plane of the
                    // face.  This is cooked (illogical!) for fastest
                    // surface walk.
                    //
                    point_type missPlanePoint =
                        curPoint + curProjectionDir*curHit.distance();

                    const labelList& masterNbrs = masterFaceFaces[curFace];

                    sqrDistance =
                        magSqr(missPlanePoint - masterFaceCentres[curFace]);

                    forAll(masterNbrs, nbrI)
                    {
                        if
                        (
                            magSqr
                            (
                                missPlanePoint
                              - masterFaceCentres[masterNbrs[nbrI]]
                            )
                         <= sqrDistance
                        )
                        {
                            closer = true;
                            curFace = masterNbrs[nbrI];
                        }
                    }

                    if (visitedTargetFace[curFace])
                    {
                        // This face has already been visited.
                        // Execute n-squared search
                        doNSquaredSearch = true;
                        break;
                    }
                }

                DebugInfo << '.';
            } while (closer);
        }

        if
        (
            doNSquaredSearch || !foundEligible
        )
        {
            nNSquaredSearches++;

            DebugInfo << "p " << curLocalPointLabel << ": ";

            result[curLocalPointLabel] = objectHit(false, -1);
            scalar minDistance = GREAT;

            forAll(masterFaces, facei)
            {
                PointHit<point_type> curHit =
                    masterFaces[facei].ray
                    (
                        curPoint,
                        curProjectionDir,
                        masterPoints,
                        alg,
                        dir
                    );

                if (curHit.hit())
                {
                    result[curLocalPointLabel] = objectHit(true, facei);
                    curFace = facei;

                    break;
                }
                else if (curHit.eligibleMiss())
                {
                    // Calculate min distance
                    scalar missDist = curHit.point().dist(curPoint);

                    if (missDist < minDistance)
                    {
                        minDistance = missDist;

                        result[curLocalPointLabel] = objectHit(false, facei);
                        curFace = facei;
                    }
                }
            }

            DebugInfo << result[curLocalPointLabel] << nl;
        }
        else
        {
            DebugInfo << 'x';
        }
    }

    DebugInfo
        << nl << "Executed " << nNSquaredSearches
        << " n-squared searches out of total of "
        << nPoints() << endl;

    return result;
}


template<class FaceList, class PointField>
template<class ToPatch>
Foam::List<Foam::objectHit>
Foam::PrimitivePatch<FaceList, PointField>::projectFaceCentres
(
    const ToPatch& targetPatch,
    const Field
    <
        typename Foam::PrimitivePatch<FaceList, PointField>::point_type
    >& projectionDirection,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // The current patch is slave, i.e. it is being projected onto the target

    if (projectionDirection.size() != this->size())
    {
        FatalErrorInFunction
            << "Projection direction field does not correspond to patch faces."
            << endl << "Size: " << projectionDirection.size()
            << " Number of points: " << this->size()
            << abort(FatalError);
    }

    labelList slaveFaceOrder = meshTools::bandCompression(faceFaces());

    // calculate master face centres
    Field<point_type> masterFaceCentres(targetPatch.size());

    const labelListList& masterFaceFaces = targetPatch.faceFaces();

    const ToPatch& masterFaces = targetPatch;

    const typename ToPatch::PointFieldType& masterPoints = targetPatch.points();

    forAll(masterFaceCentres, facei)
    {
        masterFaceCentres[facei] =
            masterFaces[facei].centre(masterPoints);
    }

    // Result
    List<objectHit> result(this->size());

    const PrimitivePatch<FaceList, PointField>& slaveFaces = *this;

    const PointField& slaveGlobalPoints = points();

    // Algorithm:
    // Loop through all points of the slave side. For every point find the
    // radius for the current contact face. If the contact point falls inside
    // the face and the radius is smaller than for all neighbouring faces,
    // the contact is found. If not, visit the neighbour closest to the
    // calculated contact point. If a single master face is visited more than
    // twice, initiate n-squared search.

    label curFace = 0;
    label nNSquaredSearches = 0;

    forAll(slaveFaceOrder, facei)
    {
        // pick up slave point and direction
        const label curLocalFaceLabel = slaveFaceOrder[facei];

        const point& curFaceCentre =
            slaveFaces[curLocalFaceLabel].centre(slaveGlobalPoints);

        const vector& curProjectionDir =
            projectionDirection[curLocalFaceLabel];

        bool closer;

        boolList visitedTargetFace(targetPatch.size(), false);
        bool doNSquaredSearch = false;

        bool foundEligible = false;

        scalar sqrDistance = GREAT;

        // Force the full search for the first point to ensure good
        // starting face
        if (facei == 0)
        {
            doNSquaredSearch = true;
        }
        else
        {
            do
            {
                closer = false;
                doNSquaredSearch = false;

                // Calculate intersection with curFace
                PointHit<point_type> curHit =
                    masterFaces[curFace].ray
                    (
                        curFaceCentre,
                        curProjectionDir,
                        masterPoints,
                        alg,
                        dir
                    );

                visitedTargetFace[curFace] = true;

                if (curHit.hit())
                {
                    result[curLocalFaceLabel] = objectHit(true, curFace);

                    break;
                }
                else
                {
                    // If a new miss is eligible, it is closer than
                    // any previous eligible miss (due to surface walk)

                    // Only grab the miss if it is eligible
                    if (curHit.eligibleMiss())
                    {
                        foundEligible = true;
                        result[curLocalFaceLabel] = objectHit(false, curFace);
                    }

                    // Find the next likely face for intersection

                    // Calculate the miss point.  This is
                    // cooked (illogical!) for fastest surface walk.
                    //
                    point_type missPlanePoint =
                        curFaceCentre + curProjectionDir*curHit.distance();

                    sqrDistance =
                        magSqr(missPlanePoint - masterFaceCentres[curFace]);

                    const labelList& masterNbrs = masterFaceFaces[curFace];

                    forAll(masterNbrs, nbrI)
                    {
                        if
                        (
                            magSqr
                            (
                                missPlanePoint
                              - masterFaceCentres[masterNbrs[nbrI]]
                            )
                         <= sqrDistance
                        )
                        {
                            closer = true;
                            curFace = masterNbrs[nbrI];
                        }
                    }

                    if (visitedTargetFace[curFace])
                    {
                        // This face has already been visited.
                        // Execute n-squared search
                        doNSquaredSearch = true;
                        break;
                    }
                }

                DebugInfo << '.';
            } while (closer);
        }

        if (doNSquaredSearch || !foundEligible)
        {
            nNSquaredSearches++;

            DebugInfo << "p " << curLocalFaceLabel << ": ";

            result[curLocalFaceLabel] = objectHit(false, -1);
            scalar minDistance = GREAT;

            forAll(masterFaces, facei)
            {
                PointHit<point_type> curHit =
                    masterFaces[facei].ray
                    (
                        curFaceCentre,
                        curProjectionDir,
                        masterPoints,
                        alg,
                        dir
                    );

                if (curHit.hit())
                {
                    result[curLocalFaceLabel] = objectHit(true, facei);
                    curFace = facei;

                    break;
                }
                else if (curHit.eligibleMiss())
                {
                    // Calculate min distance
                    scalar missDist = curHit.point().dist(curFaceCentre);

                    if (missDist < minDistance)
                    {
                        minDistance = missDist;

                        result[curLocalFaceLabel] = objectHit(false, facei);
                        curFace = facei;
                    }
                }
            }

            DebugInfo << result[curLocalFaceLabel] << nl;
        }
        else
        {
            DebugInfo << 'x';
        }
    }

    DebugInfo
        << nl
        << "Executed " << nNSquaredSearches
        << " n-squared searches out of total of "
        << this->size() << endl;

    return result;
}


// ************************************************************************* //
