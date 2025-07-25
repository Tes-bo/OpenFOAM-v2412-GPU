/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "laplacianConstraintPointSmoother.H"
#include "addToRunTimeSelectionTable.H"
#include "meshPointPatch.H"
#include "processorPointPatch.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointSmoothers
{
    defineTypeNameAndDebug(laplacianConstraintPointSmoother, 0);
    addToRunTimeSelectionTable
    (
        pointSmoother,
        laplacianConstraintPointSmoother,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointSmoothers::laplacianConstraintPointSmoother::
laplacianConstraintPointSmoother
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointSmoother(mesh, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointSmoothers::laplacianConstraintPointSmoother::calculate
(
    const labelList& facesToMove,
    const pointField& oldPoints,
    const pointField& currentPoints,
    const pointField& faceCentres,
    const vectorField& faceAreas,
    const pointField& cellCentres,
    const scalarField& cellVolumes,
    vectorField& pointDisplacement
) const
{
    // Get pointMesh so we can get at the constraints
    const auto& pMesh = pointMesh::New(mesh());

    // Number of points used in each average
    labelField counts(mesh().nPoints(), 0);

    // Reset the displacements which are about to be calculated
    reset(facesToMove, counts, pointDisplacement);

    // Get affected points
    const bitSet isMovingPoint(pointsToMove(facesToMove, true));


    // Set constraints:
    // - internal points : 0
    // - normal boundary points : 1
    // - meshPointPatch : 1 (surface)
    //                    2 (feat-edge)
    //                    3 (feat-point)
    labelList nConstraints(mesh().nPoints(), Zero);

    for (const auto& pp : pMesh.boundary())
    {
        const auto& mp = pp.meshPoints();
        const auto* mppPtr = isA<meshPointPatch>(pp);
        if (mppPtr)
        {
            const auto& pc = mppPtr->constraints();
            forAll(mp, i)
            {
                nConstraints[mp[i]] = pc[i].first();
            }
        }
        else //if (!isA<processorPointPatch>(pp)) // what about cyclic? AMI?
        {
            forAll(mp, i)
            {
                // Indirectly detect any constraint
                pointConstraint pc;
                pp.applyConstraint(i, pc);
                nConstraints[mp[i]] = pc.first();
            }
        }
    }


    // Average from equally constrained points

    const auto& edges = mesh().edges();
    const auto& pointEdges = mesh().pointEdges();
    forAll(pointEdges, pointi)
    {
        if (isMovingPoint[pointi])
        {
            const auto& pEdges = pointEdges[pointi];
            for (const label edgei : pEdges)
            {
                const label otherPointi = edges[edgei].otherVertex(pointi);
                if (nConstraints[otherPointi] >= nConstraints[pointi])
                {
                    pointDisplacement[pointi] +=
                        currentPoints[otherPointi]
                      - oldPoints[pointi];

                    ++ counts[pointi];
                }
            }
        }
    }

    // Average
    average(facesToMove, counts, pointDisplacement);


    // Make sure to set any unconnected points (or boundary feature points)
    // since otherwise they never see the effect of the boundary conditions
    forAll(counts, pointi)
    {
        if (isMovingPoint[pointi] && !counts[pointi])
        {
            pointDisplacement[pointi] = currentPoints[pointi]-oldPoints[pointi];
        }
    }
}


// ************************************************************************* //
