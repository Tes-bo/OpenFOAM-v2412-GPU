/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "pointToPointPlanarInterpolation.H"
#include "boundBox.H"
#include "Random.H"
#include "vector2D.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "OBJstream.H"
#include "Time.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToPointPlanarInterpolation, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::coordSystem::cartesian
Foam::pointToPointPlanarInterpolation::calcCoordinateSystem
(
    const pointField& points
)
{
    if (points.size() < 3)
    {
        FatalErrorInFunction
            << "Need at least 3 non-collinear points for planar interpolation,"
            << " but only had " << points.size() << " points" << nl
            << exit(FatalError);
    }

    const point& p0 = points[0];

    // Find furthest away point
    label index1 = -1;
    scalar maxDistSqr = ROOTVSMALL;

    for (label i = 1; i < points.size(); ++i)
    {
        const scalar mag2 = magSqr(points[i] - p0);

        if (maxDistSqr < mag2)
        {
            maxDistSqr = mag2;
            index1 = i;
        }
    }
    if (index1 == -1)
    {
        FatalErrorInFunction
            << "Cannot find any point that is different from first point"
            << p0 << ". Are all your points coincident?"
            << exit(FatalError);
    }

    const vector e1(normalised(points[index1] - p0));

    // Find point that is furthest perpendicular distance from the p0-p1 line
    label index2 = -1;
    maxDistSqr = ROOTVSMALL;
    for (label i = 1; i < points.size(); i++)
    {
        if (i != index1)
        {
            vector e2(points[i] - p0);
            e2.removeCollinear(e1);

            const scalar mag2 = magSqr(e2);

            if (maxDistSqr < mag2)
            {
                maxDistSqr = mag2;
                index2 = i;
            }
        }
    }
    if (index2 == -1)
    {
        FatalErrorInFunction
            << "Cannot find points that define a plane with a valid normal."
            << nl << "Have so far points " << p0 << " and " << points[index1]
            << ". Are all your points on a single line instead of a plane?"
            << exit(FatalError);
    }

    const vector n = normalised(e1 ^ (points[index2]-p0));

    DebugInFunction
        << " Used points "
        << p0 << ' ' << points[index1] << ' ' << points[index2]
        << " to define coordinate system with normal " << n << endl;

    return coordSystem::cartesian
    (
        p0,  // origin
        n,   // normal
        e1   // 0-axis
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointToPointPlanarInterpolation::calcWeights
(
    const pointField& sourcePoints,
    const pointField& destPoints
)
{
    if (nearestOnly_)
    {
        labelList destToSource;
        bool fullMatch = matchPoints
        (
            destPoints,
            sourcePoints,
            scalarField(destPoints.size(), GREAT),
            true,       // verbose
            destToSource
        );

        if (!fullMatch)
        {
            FatalErrorInFunction
                << "Did not find a corresponding sourcePoint for every face"
                << " centre" << exit(FatalError);
        }

        nearestVertex_.resize(destPoints.size());
        nearestVertexWeight_.resize(destPoints.size());
        forAll(nearestVertex_, i)
        {
            nearestVertex_[i][0] = destToSource[i];
            nearestVertex_[i][1] = -1;
            nearestVertex_[i][2] = -1;

            nearestVertexWeight_[i][0] = 1.0;
            nearestVertexWeight_[i][1] = 0.0;
            nearestVertexWeight_[i][2] = 0.0;
        }

        if (debug)
        {
            forAll(destPoints, i)
            {
                label v0 = nearestVertex_[i][0];

                Pout<< "For location " << destPoints[i]
                    << " sampling vertex " << v0
                    << " at:" << sourcePoints[v0]
                    << " distance:" << mag(sourcePoints[v0]-destPoints[i])
                    << endl;
            }

            OBJstream str
            (
                "destToSource_" + Foam::name(UPstream::myProcNo()) + ".obj"
            );
            Pout<< "pointToPointPlanarInterpolation::calcWeights :"
                << " Dumping lines from face centres to original points to "
                << str.name() << endl;

            forAll(destPoints, i)
            {
                const label v0 = nearestVertex_[i][0];
                str.writeLine(destPoints[i], sourcePoints[v0]);
            }
        }
    }
    else
    {
        auto tlocalVertices = referenceCS_.localPosition(sourcePoints);
        auto& localVertices = tlocalVertices.ref();

        const boundBox bb(localVertices, true);
        const point bbMid(bb.centre());

        DebugInFunction
            << " Perturbing points with " << perturb_
            << " fraction of a random position inside " << bb
            << " to break any ties on regular meshes." << nl
            << endl;

        Random rndGen(123456);
        forAll(localVertices, i)
        {
            localVertices[i] +=
                perturb_
               *(rndGen.position(bb.min(), bb.max())-bbMid);
        }

        // Determine triangulation
        List<vector2D> localVertices2D(localVertices.size());
        forAll(localVertices, i)
        {
            localVertices2D[i][0] = localVertices[i][0];
            localVertices2D[i][1] = localVertices[i][1];
        }

        triSurface s(triSurfaceTools::delaunay2D(localVertices2D));

        auto tlocalFaceCentres = referenceCS_.localPosition(destPoints);
        const pointField& localFaceCentres = tlocalFaceCentres();

        if (debug)
        {
            fileName outName
            (
                "triangulation_" + Foam::name(UPstream::myProcNo()) + ".obj"
            );

            Pout<< "pointToPointPlanarInterpolation::calcWeights :"
                <<" Dumping triangulated surface to " << outName << endl;

            s.write(outName);

            OBJstream os
            (
                "localFaceCentres_" + Foam::name(UPstream::myProcNo()) + ".obj"
            );
            Pout<< "pointToPointPlanarInterpolation::calcWeights :"
                << " Dumping face centres to " << os.name() << endl;

            os.write(localFaceCentres);
        }

        // Determine interpolation onto face centres.
        triSurfaceTools::calcInterpolationWeights
        (
            s,
            localFaceCentres,   // points to interpolate to
            nearestVertex_,
            nearestVertexWeight_
        );

        if (debug)
        {
            forAll(sourcePoints, i)
            {
                Pout<< "source:" << i << " at:" << sourcePoints[i]
                    << " 2d:" << localVertices[i]
                    << endl;
            }

            OBJstream str
            (
                "stencil_" + Foam::name(UPstream::myProcNo()) + ".obj"
            );
            Pout<< "pointToPointPlanarInterpolation::calcWeights :"
                << " Dumping stencil to " << str.name() << endl;

            forAll(destPoints, i)
            {
                const label v0 = nearestVertex_[i][0];
                const label v1 = nearestVertex_[i][1];
                const label v2 = nearestVertex_[i][2];

                Pout<< "For location " << destPoints[i]
                    << " 2d:" << localFaceCentres[i]
                    << " sampling vertices" << nl
                    << "    " << v0
                    << " at:" << sourcePoints[v0]
                    << " weight:" << nearestVertexWeight_[i][0] << nl;

                str.writeLine(destPoints[i], sourcePoints[v0]);

                if (v1 != -1)
                {
                    Pout<< "    " << v1
                        << " at:" << sourcePoints[v1]
                        << " weight:" << nearestVertexWeight_[i][1] << nl;
                    str.writeLine(destPoints[i], sourcePoints[v1]);
                }
                if (v2 != -1)
                {
                    Pout<< "    " << v2
                        << " at:" << sourcePoints[v2]
                        << " weight:" << nearestVertexWeight_[i][2] << nl;
                    str.writeLine(destPoints[i], sourcePoints[v2]);
                }

                Pout<< endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointToPointPlanarInterpolation::pointToPointPlanarInterpolation
(
    const pointField& sourcePoints,
    const pointField& destPoints,
    const scalar perturb,
    const bool nearestOnly
)
:
    perturb_(perturb),
    nearestOnly_(nearestOnly),
    referenceCS_(),
    nPoints_(sourcePoints.size())
{
    if (!nearestOnly)
    {
        referenceCS_ = calcCoordinateSystem(sourcePoints);
    }
    calcWeights(sourcePoints, destPoints);
}


Foam::pointToPointPlanarInterpolation::pointToPointPlanarInterpolation
(
    const coordinateSystem& referenceCS,
    const pointField& sourcePoints,
    const pointField& destPoints,
    const scalar perturb
)
:
    perturb_(perturb),
    nearestOnly_(false),
    referenceCS_(referenceCS),
    nPoints_(sourcePoints.size())
{
    calcWeights(sourcePoints, destPoints);
}


Foam::pointToPointPlanarInterpolation::pointToPointPlanarInterpolation
(
    const scalar perturb,
    const bool nearestOnly,
    const coordinateSystem& referenceCS,
    const label sourceSize,
    List<FixedList<label, 3>>&& nearestVertex,
    List<FixedList<scalar, 3>>&& nearestVertexWeight
)
:
    perturb_(perturb),
    nearestOnly_(nearestOnly),
    referenceCS_(referenceCS),
    nPoints_(sourceSize),
    nearestVertex_(std::move(nearestVertex)),
    nearestVertexWeight_(std::move(nearestVertexWeight))
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::wordList Foam::pointToPointPlanarInterpolation::timeNames
(
    const instantList& times
)
{
    wordList names(times.size());

    forAll(times, i)
    {
        names[i] = times[i].name();
    }
    return names;
}


// ************************************************************************* //
