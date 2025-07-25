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

#include "displacementMethodvolumetricBSplinesMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(displacementMethodvolumetricBSplinesMotionSolver, 0);
addToRunTimeSelectionTable
(
    displacementMethod,
    displacementMethodvolumetricBSplinesMotionSolver,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacementMethodvolumetricBSplinesMotionSolver::
displacementMethodvolumetricBSplinesMotionSolver
(
    fvMesh& mesh,
    const labelList& patchIDs
)
:
    displacementMethod(mesh, patchIDs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void displacementMethodvolumetricBSplinesMotionSolver::setMotionField
(
    const pointVectorField& pointMovement
)
{
    NotImplemented;
}


void displacementMethodvolumetricBSplinesMotionSolver::setMotionField
(
    const volVectorField& cellMovement
)
{
    NotImplemented;
}


void displacementMethodvolumetricBSplinesMotionSolver::setControlField
(
    const vectorField& controlField
)
{
    // refCast motionSolver to volumetricBSplinesMotionSolver in order to pass
    // control points movement. Safe since the motionSolver can only be an
    // volumetricBSplinesMotionSolver
    refCast<volumetricBSplinesMotionSolver>
        (motionPtr_()).setControlPointsMovement(controlField);
}


void displacementMethodvolumetricBSplinesMotionSolver::setControlField
(
    const scalarField& controlField
)
{
    NotImplemented;
}


void displacementMethodvolumetricBSplinesMotionSolver::boundControlField
(
    vectorField& controlField
)
{
    // refCast motionSolver to volumetricBSplinesMotionSolver in order to bound
    // control points movement. Safe since the motionSolver can only be an
    // volumetricBSplinesMotionSolver
    refCast<volumetricBSplinesMotionSolver>
        (motionPtr_()).boundControlPointMovement(controlField);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
