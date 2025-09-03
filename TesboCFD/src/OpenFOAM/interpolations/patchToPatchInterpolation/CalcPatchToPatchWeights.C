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

#include "PatchToPatchInterpolation.H"
#include "objectHit.H"
#include "pointHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
void PatchToPatchInterpolation<FromPatch, ToPatch>::calcPointAddressing() const
{
    // Calculate pointWeights

    pointWeightsPtr_.reset(new FieldField<Field, scalar>(toPatch_.nPoints()));
    auto& pointWeights = *pointWeightsPtr_;

    pointDistancePtr_.reset(new scalarField(toPatch_.nPoints(), GREAT));
    auto& pointDistance = *pointDistancePtr_;

    const pointField& fromPatchPoints = fromPatch_.localPoints();
    const auto& fromPatchFaces = fromPatch_.localFaces();

    const pointField& toPatchPoints = toPatch_.localPoints();
    const vectorField& projectionDirection = toPatch_.pointNormals();
    const edgeList& toPatchEdges = toPatch_.edges();
    const labelListList& toPatchPointEdges = toPatch_.pointEdges();

    if (debug)
    {
        Info<< "projecting points" << endl;
    }

    List<objectHit> proj =
        toPatch_.projectPoints(fromPatch_, projectionDirection, alg_, dir_);

    pointAddressingPtr_.reset(new labelList(proj.size(), -1));
    auto& pointAddressing = *pointAddressingPtr_;

    bool doWeights = false;

    forAll(pointAddressing, pointi)
    {
        doWeights = false;

        const auto& hitFace = fromPatchFaces[proj[pointi].hitObject()];

        point hitPoint = Zero;

        if (proj[pointi].hit())
        {
            // A hit exists
            doWeights = true;

            pointAddressing[pointi] = proj[pointi].hitObject();

            pointHit curHit =
                hitFace.ray
                (
                    toPatchPoints[pointi],
                    projectionDirection[pointi],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            // Grab distance to target
            if (dir_ == intersection::CONTACT_SPHERE)
            {
                pointDistance[pointi] =
                    hitFace.contactSphereDiameter
                    (
                        toPatchPoints[pointi],
                        projectionDirection[pointi],
                        fromPatchPoints
                    );
            }
            else
            {
                pointDistance[pointi] = curHit.distance();
            }

            // Grab hit point
            hitPoint = curHit.hitPoint();
        }
        else if (projectionTol_ > SMALL)
        {
            // Check for a near miss
            pointHit ph =
                hitFace.ray
                (
                    toPatchPoints[pointi],
                    projectionDirection[pointi],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            scalar dist =
                Foam::mag
                (
                    toPatchPoints[pointi]
                  + projectionDirection[pointi]*ph.distance()
                  - ph.missPoint()
                );

            // Calculate the local tolerance
            scalar minEdgeLength = GREAT;

            // Do shortest edge of hit object
            edgeList hitFaceEdges =
                fromPatchFaces[proj[pointi].hitObject()].edges();

            forAll(hitFaceEdges, edgeI)
            {
                minEdgeLength =
                    min
                    (
                        minEdgeLength,
                        hitFaceEdges[edgeI].mag(fromPatchPoints)
                    );
            }

            const labelList& curEdges = toPatchPointEdges[pointi];

            forAll(curEdges, edgeI)
            {
                minEdgeLength =
                    min
                    (
                        minEdgeLength,
                        toPatchEdges[curEdges[edgeI]].mag(toPatchPoints)
                    );
            }

            if (dist < minEdgeLength*projectionTol_)
            {
                // This point is being corrected
                doWeights = true;

                pointAddressing[pointi] = proj[pointi].hitObject();

                // Grab nearest point on face as hit point
                hitPoint = ph.missPoint();

                // Grab distance to target
                if (dir_ == intersection::CONTACT_SPHERE)
                {
                    pointDistance[pointi] =
                        hitFace.contactSphereDiameter
                        (
                            toPatchPoints[pointi],
                            projectionDirection[pointi],
                            fromPatchPoints
                        );
                }
                else
                {
                    pointDistance[pointi] =
                        (
                            projectionDirection[pointi]
                            /mag(projectionDirection[pointi])
                        )
                      & (hitPoint - toPatchPoints[pointi]);
                }
            }
        }

        if (doWeights)
        {
            // Set interpolation pointWeights
            const pointField hitFacePoints
            (
                hitFace.points(fromPatchPoints)
            );

            auto& pointiWeights =
                pointWeights.emplace_set(pointi, hitFacePoints.size());

            scalar sumWeight = 0;

            forAll(hitFacePoints, masterPointi)
            {
                const point& p = hitFacePoints[masterPointi];

                const scalar w =
                (
                    1.0 / (hitPoint.dist(p) + VSMALL)
                );

                pointiWeights[masterPointi] = w;
                sumWeight += w;
            }

            pointiWeights /= sumWeight;
        }
        else
        {
            pointWeights.emplace_set(pointi);
        }
    }
}


template<class FromPatch, class ToPatch>
void PatchToPatchInterpolation<FromPatch, ToPatch>::calcFaceAddressing() const
{
    faceWeightsPtr_.reset(new FieldField<Field, scalar>(toPatch_.size()));
    auto& faceWeights = *faceWeightsPtr_;

    faceDistancePtr_.reset(new scalarField(toPatch_.size(), GREAT));
    auto& faceDistance = *faceDistancePtr_;

    if (debug)
    {
        Info<< "projecting face centres" << endl;
    }

    const pointField& fromPatchPoints = fromPatch_.points();
    const auto& fromPatchFaces = fromPatch_;
    const labelListList& fromPatchFaceFaces = fromPatch_.faceFaces();

    vectorField fromPatchFaceCentres(fromPatchFaces.size());

    forAll(fromPatchFaceCentres, facei)
    {
        fromPatchFaceCentres[facei] =
            fromPatchFaces[facei].centre(fromPatchPoints);
    }

    const pointField& toPatchPoints = toPatch_.points();
    const auto& toPatchFaces = toPatch_;

    const vectorField& projectionDirection = toPatch_.faceNormals();

    List<objectHit> proj =
        toPatch_.projectFaceCentres
        (
            fromPatch_,
            projectionDirection,
            alg_,
            dir_
        );

    faceAddressingPtr_.reset(new labelList(proj.size(), -1));
    auto& faceAddressing = *faceAddressingPtr_;

    forAll(faceAddressing, facei)
    {
        if (proj[facei].hit())
        {
            // A hit exists
            faceAddressing[facei] = proj[facei].hitObject();

            const auto& hitFace = fromPatchFaces[faceAddressing[facei]];

            pointHit curHit =
                hitFace.ray
                (
                    toPatchFaces[facei].centre(toPatchPoints),
                    projectionDirection[facei],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            // grab distance to target
            faceDistance[facei] = curHit.distance();

            // grab face centre of the hit face
            const point& hitFaceCentre =
                fromPatchFaceCentres[faceAddressing[facei]];

            // grab neighbours of hit face
            const labelList& neighbours =
                fromPatchFaceFaces[faceAddressing[facei]];

            const point& hitPoint = curHit.hitPoint();

            scalar m = hitPoint.dist(hitFaceCentre);

            if
            (
                m < directHitTol                            // Direct hit
             || neighbours.empty()
            )
            {
                auto& faceiWeights = faceWeights.emplace_set(facei, 1);
                faceiWeights[0] = 1.0;
            }
            else
            {
                // set interpolation faceWeights

                // The first coefficient corresponds to the centre face.
                // The rest is ordered in the same way as the faceFaces list.

                auto& faceiWeights =
                    faceWeights.emplace_set(facei, neighbours.size() + 1);

                faceiWeights[0] = 1.0/m;

                scalar sumWeight = faceiWeights[0];

                forAll(neighbours, nbri)
                {
                    const point& p = fromPatchFaceCentres[neighbours[nbri]];

                    const scalar w =
                    (
                        1.0 / (hitPoint.dist(p) + VSMALL)
                    );

                    faceWeights[nbri+1] = w;
                    sumWeight += w;
                }

                faceiWeights /= sumWeight;
            }
        }
        else
        {
            faceWeights.emplace_set(facei);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
