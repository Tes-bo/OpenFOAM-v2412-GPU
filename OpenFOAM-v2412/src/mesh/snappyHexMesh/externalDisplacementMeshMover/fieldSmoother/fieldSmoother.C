/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

#include "fieldSmoother.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldSmoother, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldSmoother::fieldSmoother(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldSmoother::~fieldSmoother()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldSmoother::smoothNormals
(
    const label nIter,
    const bitSet& isMeshMasterPoint,
    const bitSet& isMeshMasterEdge,
    const labelList& fixedPoints,
    pointVectorField& normals
) const
{
    // Get smoothly varying internal normals field.
    Info<< typeName
        << " : Smoothing normals in interior ..." << endl;

    const edgeList& edges = mesh_.edges();

    // Points that do not change.
    bitSet isFixedPoint(mesh_.nPoints());

    // Internal points that are fixed
    forAll(fixedPoints, i)
    {
        label meshPointI = fixedPoints[i];
        isFixedPoint.set(meshPointI);
    }

    // Make sure that points that are coupled to meshPoints but not on a patch
    // are fixed as well
    syncTools::syncPointList(mesh_, isFixedPoint, maxEqOp<unsigned int>(), 0);


    // Correspondence between local edges/points and mesh edges/points
    const labelList meshPoints(identity(mesh_.nPoints()));

    // Calculate inverse sum of weights

    scalarField edgeWeights(mesh_.nEdges());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh_,
        isMeshMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );

    vectorField average;
    for (label iter = 0; iter < nIter; iter++)
    {
        meshRefinement::weightedSum
        (
            mesh_,
            isMeshMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            normals,
            average
        );
        average *= invSumWeight;

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isMeshMasterPoint,
                mag(normals-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }

        // Transfer to normals vector field
        forAll(average, pointI)
        {
            if (!isFixedPoint.test(pointI))
            {
                //full smoothing neighbours + point value
                average[pointI] = 0.5*(normals[pointI]+average[pointI]);
                normals[pointI] = normalised(average[pointI]);
            }
        }
    }
}


void Foam::fieldSmoother::smoothPatchNormals
(
    const label nIter,
    const bitSet& isPatchMasterPoint,
    const bitSet& isPatchMasterEdge,
    const indirectPrimitivePatch& adaptPatch,
    pointField& normals
) const
{
    const edgeList& edges = adaptPatch.edges();
    const labelList& meshPoints = adaptPatch.meshPoints();

    // Get smoothly varying internal normals field.
    Info<< typeName << " : Smoothing normals ..." << endl;

    scalarField edgeWeights(edges.size());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh_,
        isPatchMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );


    vectorField average;
    for (label iter = 0; iter < nIter; iter++)
    {
        meshRefinement::weightedSum
        (
            mesh_,
            isPatchMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            normals,
            average
        );
        average *= invSumWeight;

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isPatchMasterPoint,
                mag(normals-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }

        // Transfer to normals vector field
        forAll(average, pointI)
        {
            // full smoothing neighbours + point value
            average[pointI] = 0.5*(normals[pointI]+average[pointI]);
            normals[pointI] = normalised(average[pointI]);
        }
    }
}


void Foam::fieldSmoother::smoothLambdaMuDisplacement
(
    const label nIter,
    const bitSet& isMeshMasterPoint,
    const bitSet& isMeshMasterEdge,
    const bitSet& isToBeSmoothed,
    vectorField& displacement
) const
{
    const edgeList& edges = mesh_.edges();

    // Correspondence between local edges/points and mesh edges/points
    const labelList meshPoints(identity(mesh_.nPoints()));

    // Calculate inverse sum of weights
    scalarField edgeWeights(mesh_.nEdges());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh_,
        isMeshMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );

    // Get smoothly varying patch field.
    Info<< typeName << " : Smoothing displacement ..." << endl;

    const scalar lambda = 0.33;
    const scalar mu = -0.34;

    vectorField average;

    for (label iter = 0; iter < nIter; iter++)
    {
        meshRefinement::weightedSum
        (
            mesh_,
            isMeshMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            displacement,
            average
        );
        average *= invSumWeight;

        forAll(displacement, i)
        {
            if (isToBeSmoothed.test(i))
            {
                displacement[i] = (1-lambda)*displacement[i]+lambda*average[i];
            }
        }

        meshRefinement::weightedSum
        (
            mesh_,
            isMeshMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            displacement,
            average
        );
        average *= invSumWeight;


        forAll(displacement, i)
        {
            if (isToBeSmoothed.test(i))
            {
                displacement[i] = (1-mu)*displacement[i]+mu*average[i];
            }
        }


        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isMeshMasterPoint,
                mag(displacement-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }
    }
}


// ************************************************************************* //
