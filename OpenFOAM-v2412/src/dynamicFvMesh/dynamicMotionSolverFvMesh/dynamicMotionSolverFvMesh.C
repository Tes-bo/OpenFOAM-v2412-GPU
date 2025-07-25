/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "dynamicMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverFvMesh::dynamicMotionSolverFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit),
    motionPtr_(nullptr)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicMotionSolverFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    motionPtr_ = motionSolver::New(*this);
    return true;
}


Foam::dynamicMotionSolverFvMesh::dynamicMotionSolverFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    dynamicFvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    motionPtr_(motionSolver::New(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::motionSolver& Foam::dynamicMotionSolverFvMesh::motion() const
{
    return *motionPtr_;
}


bool Foam::dynamicMotionSolverFvMesh::update()
{
    fvMesh::movePoints(motionPtr_->newPoints());

    volVectorField* Uptr = getObjectPtr<volVectorField>("U");

    if (Uptr)
    {
        Uptr->correctBoundaryConditions();
    }

    return true;
}


// ************************************************************************* //
