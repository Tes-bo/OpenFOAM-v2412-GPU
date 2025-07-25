/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "volumetricBSplinesMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volumetricBSplinesMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        volumetricBSplinesMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volumetricBSplinesMotionSolver::volumetricBSplinesMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    volBSplinesBase_
    (
        const_cast<volBSplinesBase&>
        (
            volBSplinesBase::New(refCast<fvMesh>(const_cast<polyMesh&>(mesh)))
        )
    ),
    controlPointsMovement_
    (
        volBSplinesBase_.getTotalControlPointsNumber(),
        Zero
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::volumetricBSplinesMotionSolver::curPoints() const
{
    auto tPointMovement = tmp<vectorField>::New(mesh().points());
    auto& pointMovement = tPointMovement.ref();

    label pastControlPoints(0);
    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
    forAll(boxes, iNURB)
    {
        const label nb = boxes[iNURB].getControlPoints().size();
        vectorField localControlPointsMovement(nb, Zero);

        forAll(localControlPointsMovement, iCP)
        {
            localControlPointsMovement[iCP] =
                controlPointsMovement_[pastControlPoints + iCP];
        }

        tmp<vectorField>
            partialMovement
            (
                boxes[iNURB].computeNewPoints
                (
                    localControlPointsMovement
                )
            );

          pointMovement += partialMovement() - mesh().points();

          pastControlPoints += nb;
    }

    return tPointMovement;
}


void Foam::volumetricBSplinesMotionSolver::solve()
{
    // Do nothing
}


void Foam::volumetricBSplinesMotionSolver::movePoints(const pointField&)
{
    // Do nothing
}


void Foam::volumetricBSplinesMotionSolver::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::volumetricBSplinesMotionSolver::setControlPointsMovement
(
    const vectorField& controlPointsMovement
)
{
    if (controlPointsMovement_.size() != controlPointsMovement.size())
    {
        FatalErrorInFunction
           << "Attempting to replace controlPointsMovement with a set of "
           << "different size"
           << exit(FatalError);
    }
    controlPointsMovement_ = controlPointsMovement;
}


void Foam::volumetricBSplinesMotionSolver::boundControlPointMovement
(
    vectorField& controlPointsMovement
)
{
    volBSplinesBase_.boundControlPointMovement(controlPointsMovement);
}


// ************************************************************************* //
