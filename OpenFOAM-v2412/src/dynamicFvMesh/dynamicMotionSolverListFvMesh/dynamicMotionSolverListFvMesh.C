/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "dynamicMotionSolverListFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "pointMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverListFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverListFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverListFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverListFvMesh::dynamicMotionSolverListFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit),
    motionSolvers_()
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicMotionSolverListFvMesh::init
(
    const bool doInit,
    const bool mandatory
)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    IOobject ioDict
    (
        "dynamicMeshDict",
        time().constant(),
        *this,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    IOdictionary dict(ioDict);

    label i = 0;

    const auto* dictptr = dict.findDict("solvers");
    if (dictptr)
    {
        const dictionary& solverDict = *dictptr;

        motionSolvers_.setSize(solverDict.size());

        for (const entry& dEntry : solverDict)
        {
            if (dEntry.isDict())
            {
                IOobject io(ioDict);
                io.readOpt(IOobject::NO_READ);
                io.writeOpt(IOobject::AUTO_WRITE);
                io.rename(dEntry.dict().dictName());

                IOdictionary IOsolverDict
                (
                    io,
                    dEntry.dict()
                );

                motionSolvers_.set
                (
                    i++,
                    motionSolver::New(*this, IOsolverDict)
                );
            }
        }
        motionSolvers_.setSize(i);
    }
    else if (mandatory)
    {
        motionSolvers_.setSize(1);
        motionSolvers_.set(i++, motionSolver::New(*this, dict));
    }

    // Assume something changed
    return true;
}


bool Foam::dynamicMotionSolverListFvMesh::init(const bool doInit)
{
    // Fall-back to always constructing motionSolver
    return init(doInit, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverListFvMesh::~dynamicMotionSolverListFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dynamicMotionSolverListFvMesh::mapFields
(
    const mapPolyMesh& mpm
)
{
    dynamicFvMesh::mapFields(mpm);

    // Update the motionSolvers for any topo change ...
    for (auto& ms : motionSolvers_)
    {
        ms.updateMesh(mpm);
    }
}


bool Foam::dynamicMotionSolverListFvMesh::update()
{
    if (motionSolvers_.size())
    {
        // Accumulated displacement
        pointField disp(motionSolvers_[0].newPoints() - fvMesh::points());

        for (label i = 1; i < motionSolvers_.size(); i++)
        {
            disp += motionSolvers_[i].newPoints() - fvMesh::points();
        }

        fvMesh::movePoints(points() + disp);

        volVectorField* Uptr = getObjectPtr<volVectorField>("U");

        if (Uptr)
        {
            Uptr->correctBoundaryConditions();
        }
    }

    return true;
}


// ************************************************************************* //
