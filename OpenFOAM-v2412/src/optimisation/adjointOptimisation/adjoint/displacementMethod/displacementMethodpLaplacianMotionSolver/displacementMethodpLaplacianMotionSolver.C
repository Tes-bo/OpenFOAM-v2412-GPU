/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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

#include "displacementMethodpLaplacianMotionSolver.H"
#include "pLaplacianMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(displacementMethodpLaplacianMotionSolver, 1);
addToRunTimeSelectionTable
(
    displacementMethod,
    displacementMethodpLaplacianMotionSolver,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacementMethodpLaplacianMotionSolver::
displacementMethodpLaplacianMotionSolver
(
    fvMesh& mesh,
    const labelList& patchIDs
)
:
    displacementMethod(mesh, patchIDs),
    pointMotionU_(refCast<pLaplacianMotionSolver>(motionPtr_()).pointMotionU()),
    cellMotionU_(refCast<pLaplacianMotionSolver>(motionPtr_()).cellMotionU()),
    resetFields_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        ).subDict("pLaplacianMotionSolverCoeffs").getOrDefault<bool>
        (
            "resetFields",
            true
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool displacementMethodpLaplacianMotionSolver::preferPointField() const
{
    return false;
}


void displacementMethodpLaplacianMotionSolver::setMotionField
(
    const pointVectorField& pointMovement
)
{
    if (resetFields_)
    {
        pointMotionU_.primitiveFieldRef() = Zero;
        cellMotionU_.primitiveFieldRef() = Zero;
        cellMotionU_.correctBoundaryConditions();
    }

    maxDisplacement_ = SMALL;

    // Set boundary mesh movement and calculate
    // max current boundary displacement
    for (label patchI : patchIDs_)
    {
        // Set boundary field. Needed for the motionSolver class
        pointMotionU_.boundaryFieldRef()[patchI] ==
            pointMovement.boundaryField()[patchI].patchInternalField()();

        // Set boundary field values as seen from the internalField!
        // Needed for determining the max displacement
        pointMotionU_.boundaryFieldRef()[patchI].setInInternalField
        (
            pointMotionU_.primitiveFieldRef(),
            pointMovement.boundaryField()[patchI].patchInternalField()()
        );

        // Find max value
        maxDisplacement_ =
            max
            (
                maxDisplacement_,
                gMax
                (
                    mag
                    (
                        pointMotionU_.boundaryField()[patchI].
                            patchInternalField()
                    )
                )
            );
    }
    // Transfer movement to cellMotionU
    refCast<pLaplacianMotionSolver>(motionPtr_()).setBoundaryConditions();
}


void displacementMethodpLaplacianMotionSolver::setMotionField
(
    const volVectorField& cellMovement
)
{
    if (resetFields_)
    {
        pointMotionU_.primitiveFieldRef() = Zero;
        cellMotionU_.primitiveFieldRef() = Zero;
        cellMotionU_.correctBoundaryConditions();
    }

    auto& cellMotionUbf = cellMotionU_.boundaryFieldRef();
    // Set boundary mesh movement and calculate max current boundary
    // displacement
    for (const label patchI : patchIDs_)
    {
        cellMotionUbf[patchI] == cellMovement.boundaryField()[patchI];

        // Find max value
        maxDisplacement_ =
            max
            (
                maxDisplacement_,
                gMax(mag(cellMotionUbf[patchI]))
            );
    }
}


void displacementMethodpLaplacianMotionSolver::setControlField
(
    const vectorField& controlField
)
{
    NotImplemented;
}


void displacementMethodpLaplacianMotionSolver::setControlField
(
    const scalarField& controlField
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
