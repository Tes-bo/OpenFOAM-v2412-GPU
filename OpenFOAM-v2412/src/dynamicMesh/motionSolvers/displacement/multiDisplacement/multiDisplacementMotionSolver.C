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

\*----------------------------------------------------------------------------*/

#include "multiDisplacementMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiDisplacementMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        multiDisplacementMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        multiDisplacementMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiDisplacementMotionSolver::multiDisplacementMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    curPoints_(mesh.points())
{
    // Make pointDisplacement is not registered since all lower levels
    // have a pointDisplacement as well.
    pointDisplacement().checkOut();

    label i = 0;

    const dictionary& solverDict = dict.subDict("solvers");

    motionSolvers_.setSize(solverDict.size());

    for (const entry& dEntry : solverDict)
    {
        if (dEntry.isDict())
        {
            IOobject io(dict);
            io.readOpt(IOobject::NO_READ);
            io.writeOpt(IOobject::AUTO_WRITE);
            io.rename(dEntry.dict().dictName());

            IOdictionary IOsolverDict
            (
                io,
                dEntry.dict()
            );

            auto* msPtr = motionSolver::New(mesh, IOsolverDict).ptr();

            motionSolvers_.set
            (
                i,
                dynamic_cast<displacementMotionSolver*>(msPtr)
            );

            // Avoid conflicts with multiple registrations
            motionSolvers_[i].pointDisplacement().checkOut();

            i++;
        }
    }
    motionSolvers_.setSize(i);

    if (i == 0)
    {
        FatalErrorInFunction << "No displacementMotionSolvers in dictionary "
            << dict << exit(FatalError);
    }

    // Re-register so only our 'pointDisplacement' is on the database.
    pointDisplacement().checkIn();
}


Foam::multiDisplacementMotionSolver::
multiDisplacementMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    curPoints_(mesh.points())
{
    // Make pointDisplacement is not registered since all lower levels
    // have a pointDisplacement as well.
    this->pointDisplacement().checkOut();

    label i = 0;

    const dictionary& solverDict = dict.subDict("solvers");

    motionSolvers_.setSize(solverDict.size());

    for (const entry& dEntry : solverDict)
    {
        if (dEntry.isDict())
        {
            IOobject io(dict);
            io.readOpt(IOobject::NO_READ);
            io.writeOpt(IOobject::AUTO_WRITE);
            io.rename(dEntry.dict().dictName());

            IOdictionary IOsolverDict
            (
                io,
                dEntry.dict()
            );

            auto msPtr = displacementMotionSolver::New
            (
                dEntry.keyword(),
                mesh,
                IOsolverDict,
                pointDisplacement,
                points0
            );

            // Avoid conflicts with multiple registrations
            msPtr->pointDisplacement().checkOut();

            motionSolvers_.set(i++, msPtr);
        }
    }
    motionSolvers_.setSize(i);

    if (i == 0)
    {
        FatalErrorInFunction << "No displacementMotionSolvers in dictionary "
            << dict << exit(FatalError);
    }

    // Re-register so only our 'pointDisplacement' is on the database.
    this->pointDisplacement().checkIn();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::multiDisplacementMotionSolver::curPoints() const
{
    return curPoints_;
}


void Foam::multiDisplacementMotionSolver::solve()
{
    if (!motionSolvers_.size())
    {
        return;
    }

    // Bit tricky:
    // - make sure only one copy of pointDisplacement is registered. This is
    //   for if cellDisplacement tries to look up the pointDisplacement it
    //   looks up its version. Or maybe we always should use our own version
    //   only?
    // - copy the last set of calculated points into our copy (curPoints_)
    // - move the mesh to update the faceCentres, cellCentres etc. This assumes
    //   that we can call movePoints() multiple times inside a time step.
    //   (note that this is supported in pimpleFoam with the
    //    moveMeshOuterCorrectors option)

    pointDisplacement().checkOut();

    // Doing first motion solver
    motionSolvers_[0].pointDisplacement().checkIn();
    // Take over my bc values
    motionSolvers_[0].pointDisplacement() == pointDisplacement();

    motionSolvers_[0].solve();
    motionSolvers_[0].pointDisplacement().checkOut();

    // Update my values
    curPoints_ = motionSolvers_[0].curPoints();
    pointDisplacement() == motionSolvers_[0].pointDisplacement();

    for (label i = 1; i < motionSolvers_.size(); i++)
    {
        // Doing other motion solvers using new locations/face/cellCentres etc.
        const_cast<polyMesh&>(mesh()).movePoints(curPoints_);
        motionSolvers_[i].pointDisplacement().checkIn();

        // Take over my bc values
        motionSolvers_[i].pointDisplacement() == pointDisplacement();

        motionSolvers_[i].solve();
        motionSolvers_[i].pointDisplacement().checkOut();

        // Update my values
        curPoints_ = motionSolvers_[i].curPoints();
        pointDisplacement() == motionSolvers_[i].pointDisplacement();
    }

    pointDisplacement().checkIn();

    // Push my pointDisplacement onto all motionSolvers
    for (auto& ms : motionSolvers_)
    {
        ms.pointDisplacement() == pointDisplacement();
    }
}


void Foam::multiDisplacementMotionSolver::movePoints
(
    const pointField& newPoints
)
{
    curPoints_ = newPoints;
    for (auto& ms : motionSolvers_)
    {
        ms.movePoints(newPoints);
    }
}


void Foam::multiDisplacementMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    for (auto& ms : motionSolvers_)
    {
        ms.updateMesh(mpm);
    }
}


// ************************************************************************* //
