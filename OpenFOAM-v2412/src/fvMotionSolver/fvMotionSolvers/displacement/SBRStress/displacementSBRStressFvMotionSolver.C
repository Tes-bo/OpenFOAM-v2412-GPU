/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "displacementSBRStressFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "mapPolyMesh.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementSBRStressFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementSBRStressFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        displacementSBRStressFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementSBRStressFvMotionSolver::displacementSBRStressFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement().dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement().boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    )
{}


Foam::displacementSBRStressFvMotionSolver::
displacementSBRStressFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            displacementMotionSolver::pointDisplacement().dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>
        (
            displacementMotionSolver::pointDisplacement().boundaryField()
        )
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementSBRStressFvMotionSolver::
~displacementSBRStressFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementSBRStressFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    // Evaluate the bcs so they are consistent with the internal field
    // Might fight the multi-patch behaviour inside volPointInterpolate
    if
    (
        pointDisplacement_.boundaryField().size()
     != cellDisplacement_.boundaryField().size()
    )
    {
        pointDisplacement_.correctBoundaryConditions();
    }

    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement().primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::displacementSBRStressFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    const surfaceScalarField Df
    (
        dimensionedScalar("viscosity", dimViscosity, 1.0)
       *diffusivityPtr_->operator()()
    );

    // Make sure the cellMotion bcs are consistent with the pointDisplacement.
    // This makes sure the grad below uses more up-to-date values.
    cellDisplacement_.boundaryFieldRef().updateCoeffs();

    volTensorField gradCd("gradCd", fvc::grad(cellDisplacement_));

    fv::options& fvOptions(fv::options::New(fvMesh_));

    fvVectorMatrix TEqn
    (
        fvm::laplacian
        (
            2*Df,
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )

      + fvc::div
        (
            Df
           *(
               fvc::dotInterpolate
               (
                   cellDisplacement_.mesh().Sf(),
                   gradCd.T() - gradCd
               )

               // Solid-body rotation "lambda" term
             - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
            )
        )

        /*
      - fvc::laplacian
        (
            2*Df,
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )

      + fvc::div
        (
            Df
           *(
               fvc::dotInterpolate
               (
                   cellDisplacement_.mesh().Sf(),
                   gradCd + gradCd.T()
               )

               // Solid-body rotation "lambda" term
             - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
           )
        )
        */
     ==
        fvOptions(cellDisplacement_)
    );

    fvOptions.constrain(TEqn);
    TEqn.solveSegregatedOrCoupled();
    fvOptions.correct(cellDisplacement_);
}


void Foam::displacementSBRStressFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


// ************************************************************************* //
