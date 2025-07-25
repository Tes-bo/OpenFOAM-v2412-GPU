/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "displacementMethodlaplacianMotionSolver.H"
#include "laplacianMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(displacementMethodlaplacianMotionSolver, 1);
addToRunTimeSelectionTable
(
    displacementMethod,
    displacementMethodlaplacianMotionSolver,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacementMethodlaplacianMotionSolver::displacementMethodlaplacianMotionSolver
(
    fvMesh& mesh,
    const labelList& patchIDs
)
:
    displacementMethod(mesh, patchIDs),
    pointMotionU_(refCast<laplacianMotionSolver>(motionPtr_()).pointMotionU()),
    cellMotionU_(refCast<laplacianMotionSolver>(motionPtr_()).cellMotionU()),
    resetFields_
    (
        IOdictionary::readContents
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            )
        ).subDict("laplacianMotionSolverCoeffs").getOrDefault<bool>
        (
            "resetFields",
            true
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool displacementMethodlaplacianMotionSolver::preferPointField() const
{
    return false;
}


void displacementMethodlaplacianMotionSolver::setMotionField
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
    for (const label patchI : patchIDs_)
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
    refCast<laplacianMotionSolver>(motionPtr_()).setBoundaryConditions();
}


void displacementMethodlaplacianMotionSolver::setMotionField
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


void displacementMethodlaplacianMotionSolver::setControlField
(
    const vectorField& controlField
)
{
    NotImplemented;
}


void displacementMethodlaplacianMotionSolver::setControlField
(
    const scalarField& controlField
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
