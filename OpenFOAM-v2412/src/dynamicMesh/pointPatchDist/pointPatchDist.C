/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "pointPatchDist.H"
#include "externalPointEdgePoint.H"
#include "pointMesh.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointPatchDist::pointPatchDist
(
    const pointMesh& pMesh,
    const labelHashSet& patchIDs,
    const pointField& points
)
:
    pointScalarField
    (
        IOobject
        (
            "pointDistance",
            pMesh.db().time().timeName(),
            pMesh.db()
        ),
        pMesh,
        dimensionedScalar("y", dimLength, GREAT)
    ),
    points_(points),
    patchIDs_(patchIDs),
    nUnset_(0)
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointPatchDist::correct()
{
    const pointBoundaryMesh& pbm = mesh().boundary();

    label nPoints = 0;

    for (const label patchi : patchIDs_)
    {
        nPoints += pbm[patchi].meshPoints().size();
    }

    externalPointEdgePoint::trackingData td(points_);

    // Set initial changed points to all the patch points(if patch present)
    List<externalPointEdgePoint> wallInfo(nPoints);
    labelList wallPoints(nPoints);
    nPoints = 0;

    for (const label patchi : patchIDs_)
    {
        // Retrieve the patch now we have its index in patches.

        const labelList& mp = pbm[patchi].meshPoints();

        forAll(mp, ppI)
        {
            label meshPointi = mp[ppI];
            wallPoints[nPoints] = meshPointi;
            wallInfo[nPoints] = externalPointEdgePoint
            (
                td.points_[meshPointi],
                0.0
            );
            nPoints++;
        }
    }

    // Current info on points
    List<externalPointEdgePoint> allPointInfo(mesh()().nPoints());

    // Current info on edges
    List<externalPointEdgePoint> allEdgeInfo(mesh()().nEdges());

    PointEdgeWave
    <
        externalPointEdgePoint,
        externalPointEdgePoint::trackingData
    > wallCalc
    (
        mesh()(),
        wallPoints,
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints(), // max iterations
        td
    );

    pointScalarField& psf = *this;


    forAll(allPointInfo, pointi)
    {
        if (allPointInfo[pointi].valid(td))
        {
            psf[pointi] = Foam::sqrt(allPointInfo[pointi].distSqr());
        }
        else
        {
            nUnset_++;
        }
    }
}


// ************************************************************************* //
