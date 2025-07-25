/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "points0MotionSolver.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(points0MotionSolver, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::points0MotionSolver::points0IO(const polyMesh& mesh)
{
    const word instance =
        mesh.time().findInstance
        (
            mesh.meshDir(),
            "points0",
            IOobject::READ_IF_PRESENT
        );
    IOobject io
    (
        "points0",
        instance,
        polyMesh::meshSubDir,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    // If points0 are located in constant directory, verify their existence
    // or fallback to a copy of the original mesh points
    if
    (
        instance == mesh.time().constant()
     && !io.typeHeaderOk<pointIOField>()
    )
    {
        io.rename("points");
    }

    return io;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::points0MotionSolver::points0MotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    zoneMotion(coeffDict(), mesh),
    points0_(points0IO(mesh))
{
    if
    (
        FieldBase::allowConstructFromLargerSize
     && (points0_.size() > mesh.nPoints())
    )
    {
        // Allowed
    }
    else if (points0_.size() != mesh.nPoints())
    {
        const fileName fName
        (
            IOobject
            (
                "points",
                time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ).typeFilePath<pointIOField>()
        );

        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file " << fName << nl
            << exit(FatalError);
    }
}


Foam::points0MotionSolver::points0MotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointIOField& points0,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    zoneMotion(coeffDict(), mesh),
    points0_(points0)
{
    if (points0_.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file " << points0.filePath()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::points0MotionSolver::movePoints(const pointField&)
{}


void Foam::points0MotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    // pointMesh already updates pointFields

    motionSolver::updateMesh(mpm);

    // Map points0_. Bit special since we somehow have to come up with
    // a sensible points0 position for introduced points.
    // Find out scaling between points0 and current points

    // Get the new points either from the map or the mesh
    const pointField& points =
    (
        mpm.hasMotionPoints()
      ? mpm.preMotionPoints()
      : mesh().points()
    );

    // Note: boundBox does reduce
    const vector span0 = boundBox(points0_).span();
    const vector span = boundBox(points).span();

    vector scaleFactors(cmptDivide(span0, span));

    pointField newPoints0(mpm.pointMap().size());

    forAll(newPoints0, pointi)
    {
        label oldPointi = mpm.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            label masterPointi = mpm.reversePointMap()[oldPointi];

            if (masterPointi == pointi)
            {
                newPoints0[pointi] = points0_[oldPointi];
            }
            else
            {
                // New point - assume motion is scaling
                newPoints0[pointi] = points0_[oldPointi] + cmptMultiply
                (
                    scaleFactors,
                    points[pointi] - points[masterPointi]
                );
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot determine coordinates of introduced vertices."
                << " New vertex " << pointi << " at coordinate "
                << points[pointi] << exit(FatalError);
        }
    }

    twoDCorrectPoints(newPoints0);

    points0_.transfer(newPoints0);

    // points0 changed - set to write and check-in to database
    points0_.rename("points0");
    points0_.writeOpt(IOobject::AUTO_WRITE);
    points0_.instance() = time().timeName();
    points0_.checkIn();
}


// ************************************************************************* //
