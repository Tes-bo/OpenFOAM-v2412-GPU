/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "pointFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointFile, 0);
    addToRunTimeSelectionTable(initialPointsMethod, pointFile, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointFile::pointFile
(
    const dictionary& initialPointsDict,
    const Time& runTime,
    Random& rndGen,
    const conformationSurfaces& geometryToConformTo,
    const cellShapeControl& cellShapeControls,
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
:
    initialPointsMethod
    (
        typeName,
        initialPointsDict,
        runTime,
        rndGen,
        geometryToConformTo,
        cellShapeControls,
        decomposition
    ),
    pointFileName_(detailsDict().get<fileName>("pointFile").expand()),
    insideOutsideCheck_(detailsDict().get<Switch>("insideOutsideCheck")),
    randomiseInitialGrid_(detailsDict().get<Switch>("randomiseInitialGrid")),
    randomPerturbationCoeff_
    (
        detailsDict().get<scalar>("randomPerturbationCoeff")
    )
{
    Info<< "    Inside/Outside check is " << insideOutsideCheck_.c_str()
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::List<Vb::Point> Foam::pointFile::initialPoints() const
{
    pointField points;
    {
        // Look for points
        IOobject pointsIO
        (
            pointFileName_.name(),
            time().timeName(),
            time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        Info<< "    Inserting points from file " << pointFileName_ << endl;

        // See if processor local file
        if (pointsIO.typeHeaderOk<pointIOField>(true))
        {
            // Found it (processor local)
            points = pointIOField(pointsIO);

            if (points.empty())
            {
                FatalErrorInFunction
                    << "Point file contain no points"
                    << exit(FatalError) << endl;
            }

            if (Pstream::parRun())
            {
                // assume that the points in each processor file
                // are unique.  They are unlikely to belong on the current
                // processor as the background mesh is unlikely to be the same.
                decomposition().distributePoints(points);
            }
        }
        else if (Pstream::parRun())
        {
            // See if points can be found in parent directory
            // (only if timeName = constant)
            points = pointIOField
            (
                IOobject
                (
                    pointFileName_.name(),
                    time().caseConstant(),
                    time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            if (points.empty())
            {
                FatalErrorInFunction
                    << "Point file contain no points"
                    << exit(FatalError) << endl;
            }

            // Points are assumed to be covering the whole
            // domain, so filter the points to be only those on this processor
            boolList procPt(decomposition().positionOnThisProcessor(points));

            List<boolList> allProcPt(Pstream::nProcs());

            allProcPt[Pstream::myProcNo()] = procPt;
            Pstream::allGatherList(allProcPt);

            forAll(procPt, ptI)
            {
                bool foundAlready = false;

                forAll(allProcPt, proci)
                {
                    // If a processor with a lower index has found this point
                    // to insert already, defer to it and don't insert.
                    if (foundAlready)
                    {
                        allProcPt[proci][ptI] = false;
                    }
                    else if (allProcPt[proci][ptI])
                    {
                        foundAlready = true;
                    }
                }
            }

            procPt = allProcPt[Pstream::myProcNo()];

            inplaceSubset(procPt, points);
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find points file " << pointsIO.objectPath()
                << exit(FatalError) << endl;
        }
    }

    Field<bool> insidePoints(points.size(), true);

    if (insideOutsideCheck_)
    {
        insidePoints = geometryToConformTo().wellInside
        (
            points,
            minimumSurfaceDistanceCoeffSqr_
           *sqr(cellShapeControls().cellSize(points))
        );
    }

    DynamicList<Vb::Point> initialPoints(insidePoints.size()/10);

    forAll(insidePoints, i)
    {
        if (insidePoints[i])
        {
            point& p = points[i];

            if (randomiseInitialGrid_)
            {
                p.x() +=
                    randomPerturbationCoeff_
                   *(rndGen().sample01<scalar>() - 0.5);
                p.y() +=
                    randomPerturbationCoeff_
                   *(rndGen().sample01<scalar>() - 0.5);
                p.z() +=
                    randomPerturbationCoeff_
                   *(rndGen().sample01<scalar>() - 0.5);
            }

            initialPoints.append(Vb::Point(p.x(), p.y(), p.z()));
        }
    }

    initialPoints.shrink();

    label nPointsRejected = points.size() - initialPoints.size();

    if (Pstream::parRun())
    {
        reduce(nPointsRejected, sumOp<label>());
    }

    if (nPointsRejected)
    {
        Info<< "    " << nPointsRejected << " points rejected from "
            << pointFileName_.name() << endl;
    }

    return initialPoints;
}


// ************************************************************************* //
