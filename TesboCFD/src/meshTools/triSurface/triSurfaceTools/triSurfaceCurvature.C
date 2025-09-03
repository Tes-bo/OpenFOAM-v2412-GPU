/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "triSurfaceTools.H"

#include "triSurface.H"
#include "MeshedSurfaces.H"
#include "triSurfaceFields.H"
#include "OFstream.H"
#include "plane.H"
#include "tensor2D.H"
#include "symmTensor2D.H"
#include "scalarMatrices.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::triSurfaceTools::vertexNormalWeight
(
    const triFace& f,
    const label pI,
    const vector& fN,
    const UList<point>& points
)
{
    label index = f.find(pI);

    if (index == -1)
    {
        FatalErrorInFunction
            << "Point not in face" << abort(FatalError);
    }

    const vector e1 = points[f[index]] - points[f[f.fcIndex(index)]];
    const vector e2 = points[f[index]] - points[f[f.rcIndex(index)]];

    return mag(fN)/(magSqr(e1)*magSqr(e2) + VSMALL);
}


Foam::tmp<Foam::vectorField>
Foam::triSurfaceTools::vertexNormals(const triSurface& surf)
{
    // Weighted average of normals of faces attached to the vertex
    // Weight = fA / (mag(e0)^2 * mag(e1)^2);

    Info<< "Calculating vertex normals" << endl;

    auto tpointNormals = tmp<vectorField>::New(surf.nPoints(), Zero);
    auto& pointNormals = tpointNormals.ref();

    const pointField& points = surf.points();
    const labelListList& pointFaces = surf.pointFaces();
    const labelList& meshPoints = surf.meshPoints();

    forAll(pointFaces, pI)
    {
        const labelList& pFaces = pointFaces[pI];

        for (const label facei : pFaces)
        {
            const triFace& f = surf[facei];

            const vector areaNorm = f.areaNormal(points);

            const scalar weight = vertexNormalWeight
            (
                f,
                meshPoints[pI],
                areaNorm,
                points
            );

            pointNormals[pI] += weight * areaNorm;
        }

        pointNormals[pI].normalise();
    }

    return tpointNormals;
}


Foam::tmp<Foam::triadField>
Foam::triSurfaceTools::vertexTriads
(
    const triSurface& surf,
    const vectorField& pointNormals
)
{
    const pointField& points = surf.points();
    const Map<label>& meshPointMap = surf.meshPointMap();

    auto tpointTriads = tmp<triadField>::New(points.size());
    auto& pointTriads = tpointTriads.ref();

    forAll(points, pI)
    {
        const point& pt = points[pI];
        const vector& normal = pointNormals[meshPointMap[pI]];

        if (mag(normal) < SMALL)
        {
            pointTriads[meshPointMap[pI]] = triad::unset;
            continue;
        }

        plane p(pt, normal);

        // Pick arbitrary point in plane
        vector dir1 = normalised(pt - p.somePointInPlane(1e-3));
        vector dir2 = normalised(dir1 ^ normal);

        pointTriads[meshPointMap[pI]] = triad(dir1, dir2, normal);
    }

    return tpointTriads;
}


Foam::tmp<Foam::scalarField>
Foam::triSurfaceTools::curvatures
(
    const triSurface& surf,
    const vectorField& pointNormals,
    const triadField& pointTriads
)
{
    Info<< "Calculating face curvature" << endl;

    const pointField& points = surf.points();
    const labelList& meshPoints = surf.meshPoints();
    const Map<label>& meshPointMap = surf.meshPointMap();

    List<symmTensor2D> pointFundamentalTensors
    (
        points.size(),
        symmTensor2D::zero
    );

    scalarList accumulatedWeights(points.size(), Zero);

    forAll(surf, fI)
    {
        const triFace& f = surf[fI];
        const edgeList fEdges = f.edges();

        // Calculate the edge vectors and the normal differences
        vectorField edgeVectors(f.size(), Zero);
        vectorField normalDifferences(f.size(), Zero);

        forAll(fEdges, feI)
        {
            const edge& e = fEdges[feI];

            edgeVectors[feI] = e.vec(points);
            normalDifferences[feI] =
               pointNormals[meshPointMap[e[0]]]
             - pointNormals[meshPointMap[e[1]]];
        }

        // Set up a local coordinate system for the face
        const vector& e0 = edgeVectors[0];
        const vector eN = f.areaNormal(points);
        const vector e1 = (e0 ^ eN);

        if (magSqr(eN) < ROOTVSMALL)
        {
            continue;
        }

        triad faceCoordSys(e0, e1, eN);
        faceCoordSys.normalize();

        // Construct the matrix to solve
        scalarSymmetricSquareMatrix T(3, 0);
        scalarDiagonalMatrix Z(3, 0);

        // Least Squares
        for (label i = 0; i < 3; ++i)
        {
            scalar x = edgeVectors[i] & faceCoordSys[0];
            scalar y = edgeVectors[i] & faceCoordSys[1];

            T(0, 0) += sqr(x);
            T(1, 0) += x*y;
            T(1, 1) += sqr(x) + sqr(y);
            T(2, 1) += x*y;
            T(2, 2) += sqr(y);

            scalar dndx = normalDifferences[i] & faceCoordSys[0];
            scalar dndy = normalDifferences[i] & faceCoordSys[1];

            Z[0] += dndx*x;
            Z[1] += dndx*y + dndy*x;
            Z[2] += dndy*y;
        }

        // Perform Cholesky decomposition and back substitution.
        // Decomposed matrix is in T and solution is in Z.
        LUsolve(T, Z);
        symmTensor2D secondFundamentalTensor(Z[0], Z[1], Z[2]);

        // Loop over the face points adding the contribution of the face
        // curvature to the points.
        forAll(f, fpI)
        {
            const label patchPointIndex = meshPointMap[f[fpI]];

            const triad& ptCoordSys = pointTriads[patchPointIndex];

            if (!ptCoordSys.set())
            {
                continue;
            }

            // Rotate faceCoordSys to ptCoordSys
            tensor rotTensor = rotationTensor(ptCoordSys[2], faceCoordSys[2]);
            triad rotatedFaceCoordSys = rotTensor & tensor(faceCoordSys);

            // Project the face curvature onto the point plane

            vector2D cmp1
            (
                ptCoordSys[0] & rotatedFaceCoordSys[0],
                ptCoordSys[0] & rotatedFaceCoordSys[1]
            );
            vector2D cmp2
            (
                ptCoordSys[1] & rotatedFaceCoordSys[0],
                ptCoordSys[1] & rotatedFaceCoordSys[1]
            );

            tensor2D projTensor
            (
                cmp1,
                cmp2
            );

            symmTensor2D projectedFundamentalTensor
            (
                projTensor.x() & (secondFundamentalTensor & projTensor.x()),
                projTensor.x() & (secondFundamentalTensor & projTensor.y()),
                projTensor.y() & (secondFundamentalTensor & projTensor.y())
            );

            // Calculate weight
            // TODO: Voronoi area weighting
            scalar weight = triSurfaceTools::vertexNormalWeight
            (
                f,
                meshPoints[patchPointIndex],
                f.areaNormal(points),
                points
            );

            // Sum contribution of face to this point
            pointFundamentalTensors[patchPointIndex] +=
                weight*projectedFundamentalTensor;

            accumulatedWeights[patchPointIndex] += weight;
        }

        if (false)
        {
            Info<< "Points = " << points[f[0]] << " "
                << points[f[1]] << " "
                << points[f[2]] << endl;
            Info<< "edgeVecs = " << edgeVectors[0] << " "
                << edgeVectors[1] << " "
                << edgeVectors[2] << endl;
            Info<< "normDiff = " << normalDifferences[0] << " "
                << normalDifferences[1] << " "
                << normalDifferences[2] << endl;
            Info<< "faceCoordSys = " << faceCoordSys << endl;
            Info<< "T = " << T << endl;
            Info<< "Z = " << Z << endl;
        }
    }

    auto tcurvatureAtPoints = tmp<scalarField>::New(points.size(), Zero);
    scalarField& curvatureAtPoints = tcurvatureAtPoints.ref();

    forAll(curvatureAtPoints, pI)
    {
        pointFundamentalTensors[pI] /= (accumulatedWeights[pI] + SMALL);

        vector2D principalCurvatures(eigenValues(pointFundamentalTensors[pI]));

        //scalar curvature =
        //    (principalCurvatures[0] + principalCurvatures[1])/2;
        scalar curvature = max
        (
            mag(principalCurvatures[0]),
            mag(principalCurvatures[1])
        );
        //scalar curvature = principalCurvatures[0]*principalCurvatures[1];

        curvatureAtPoints[meshPoints[pI]] = curvature;
    }

    return tcurvatureAtPoints;
}


Foam::tmp<Foam::scalarField>
Foam::triSurfaceTools::curvatures
(
    const triSurface& surf
)
{
    tmp<vectorField> norms = triSurfaceTools::vertexNormals(surf);
    tmp<triadField> triads = triSurfaceTools::vertexTriads(surf, norms());

    tmp<scalarField> curv = curvatures(surf, norms(), triads());
    norms.clear();
    triads.clear();

    return curv;
}


Foam::tmp<Foam::scalarField>
Foam::triSurfaceTools::writeCurvature
(
    const Time& runTime,
    const word& basename,
    const triSurface& surf
)
{
    Info<< "Extracting curvature of surface at the points." << endl;

    tmp<scalarField> tcurv = triSurfaceTools::curvatures(surf);
    scalarField& curv = tcurv.ref();

    triSurfacePointScalarField outputField
    (
        IOobject
        (
            basename + ".curvature",
            runTime.constant(),
            "triSurface",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        scalarField()
    );

    outputField.swap(curv);
    outputField.write();
    outputField.swap(curv);

    return tcurv;
}


// ************************************************************************* //
