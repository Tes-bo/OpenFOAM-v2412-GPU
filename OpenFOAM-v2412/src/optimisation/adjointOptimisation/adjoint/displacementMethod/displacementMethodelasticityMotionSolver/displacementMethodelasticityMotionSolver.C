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

#include "displacementMethodelasticityMotionSolver.H"
#include "elasticityMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(displacementMethodelasticityMotionSolver, 1);
addToRunTimeSelectionTable
(
    displacementMethod,
    displacementMethodelasticityMotionSolver,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacementMethodelasticityMotionSolver::
displacementMethodelasticityMotionSolver
(
    fvMesh& mesh,
    const labelList& patchIDs
)
:
    displacementMethod(mesh, patchIDs),
    pointMotionU_(refCast<elasticityMotionSolver>(motionPtr_()).pointMotionU()),
    cellMotionU_(refCast<elasticityMotionSolver>(motionPtr_()).cellMotionU()),
    /*
    // Getting a reference from the database is dangerous since multiple
    // fields with the same name might be registered
    pointMotionU_
    (
        const_cast<pointVectorField&>
        (
            mesh_.lookupObject<pointVectorField>("pointMotionU")
        )
    ),
    cellMotionU_
    (
        const_cast<volVectorField&>
        (
            mesh_.lookupObject<volVectorField>("cellMotionU")
        )
    )
    */
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
        ).subDict("elasticityMotionSolverCoeffs").getOrDefault<bool>
        (
            "resetFields",
            true
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool displacementMethodelasticityMotionSolver::preferPointField() const
{
    return false;
}


void displacementMethodelasticityMotionSolver::setMotionField
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

    // Update the boundary conditions of the pointField in order to make sure
    // that the boundary will move according to the initial BCs without the
    // interference of the volPointInterpolation in the elasticityMotionSolver
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

    // Update the volField boundary conditions,
    // used for the elasticity PDEs solution
    const pointField& points = mesh_.points();
    for (label patchI : patchIDs_)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        fvPatchVectorField& bField = cellMotionU_.boundaryFieldRef()[patchI];
        forAll(patch, fI)
        {
            bField[fI] = patch[fI].average(points, pointMovement);
        }
    }
}


void displacementMethodelasticityMotionSolver::setMotionField
(
    const volVectorField& cellMovement
)
{
    auto& cellMotionUbf = cellMotionU_.boundaryFieldRef();

    // Set boundary mesh movement and calculate
    // max current boundary displacement
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


void displacementMethodelasticityMotionSolver::setControlField
(
    const vectorField& controlField
)
{
    NotImplemented;
}


void displacementMethodelasticityMotionSolver::setControlField
(
    const scalarField& controlField
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
