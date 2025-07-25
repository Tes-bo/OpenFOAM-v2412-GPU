/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd
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

#include "addToRunTimeSelectionTable.H"
#include "dynamicMotionSolverTopoFvMesh.H"
#include "mapPolyMesh.H"
#include "OBJstream.H"
#include "Time.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverTopoFvMesh, 0);

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverTopoFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverTopoFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverTopoFvMesh::dynamicMotionSolverTopoFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    topoChangerFvMesh(io, doInit)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicMotionSolverTopoFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        topoChangerFvMesh::init(doInit);
    }

    motionPtr_ = motionSolver::New(*this);

    // Assume something changed
    return true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverTopoFvMesh::~dynamicMotionSolverTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMotionSolverTopoFvMesh::update()
{
    // Clear moving flag. This is currently required since geometry calculation
    // might get triggered when doing processor patches.
    // (TBD: should be in changeMesh if no inflation?)
    moving(false);
    // Do mesh changes (not using inflation - points added directly into mesh)
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(false);

    if (topoChangeMap)
    {
        Info << "Executing mesh topology update" << endl;
        motionPtr_->updateMesh(topoChangeMap());

        setV0() = V();

        pointField newPoints(motionPtr_->newPoints());
        movePoints(newPoints);

        if (debug)
        {
            OBJstream osOld("oldPts_" + time().timeName() + ".obj");
            osOld.write(oldPoints());

            OBJstream osNew("newPts_" + time().timeName() + ".obj");
            osNew.write(points());
        }
    }
    else
    {
        // Calculate the new point positions using the motion solver
        pointField newPoints(motionPtr_->newPoints());

        // The mesh now contains the cells with zero volume
        Info << "Executing mesh motion" << endl;
        movePoints(newPoints);
    }

    return true;
}


// ************************************************************************* //
