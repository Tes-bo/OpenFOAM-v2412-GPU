/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "elasticityMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "wallDist.H"
#include "fixedValuePointPatchFields.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "fvmDiv.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elasticityMotionSolver, 1);

    addToRunTimeSelectionTable
    (
        motionSolver,
        elasticityMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::elasticityMotionSolver::setBoundaryConditions()
{
    // Adjust boundary conditions based on the steps to be executed
    forAll(cellMotionU_.boundaryField(), patchI)
    {
        fvPatchVectorField& bc =
            cellMotionU_.boundaryFieldRef()[patchI];
        if (isA<fixedValueFvPatchVectorField>(bc))
        {
            auto& fixedValueBCs =
                refCast<fixedValueFvPatchVectorField>(bc);
            fixedValueBCs == fixedValueBCs/scalar(nSteps_);
        }
    }
    /*
    // Adjust boundary conditions based on the steps to be executed
    forAll(pointMotionU_.boundaryField(), patchI)
    {
        pointPatchVectorField& pointBCs =
            pointMotionU_.boundaryFieldRef()[patchI];
        if (isA<fixedValuePointPatchVectorField>(pointBCs))
        {
            auto& fixedValueBCs =
                refCast<fixedValuePointPatchVectorField>(pointBCs);
            fixedValueBCs == fixedValueBCs/scalar(nSteps_);
        }
    }

    // Copy boundary conditions to internalField
    // Needed for interpolation to faces
    pointMotionU_.boundaryFieldRef().updateCoeffs();

    // Interpolate boundary conditions from points to faces
    forAll(cellMotionU_.boundaryField(), pI)
    {
        fvPatchVectorField& bField = cellMotionU_.boundaryFieldRef()[pI];
        if (isA<fixedValueFvPatchVectorField>(bField))
        {
            const pointField& points = fvMesh_.points();
            const polyPatch& patch = mesh().boundaryMesh()[pI];
            forAll(bField, fI)
            {
                bField[fI] = patch[fI].average(points, pointMotionU_);
            }
        }
    }
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::elasticityMotionSolver::elasticityMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    fvMesh_
    (
        const_cast<fvMesh&>
        (
            refCast<const fvMesh>(mesh)
        )
    ),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedVector(dimless, Zero),
        fixedValuePointPatchVectorField::typeName
    ),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointMotionU_.dimensions(), Zero),
        pointMotionU_.boundaryField().types()
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
    ),
    nSteps_(this->coeffDict().get<label>("steps")),
    nIters_(this->coeffDict().get<label>("iters")),
    tolerance_(this->coeffDict().get<scalar>("tolerance"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::elasticityMotionSolver::curPoints() const
{
    return tmp<pointField>::New(mesh().points());
}


void Foam::elasticityMotionSolver::solve()
{
    // Re-init to zero
    cellMotionU_.primitiveFieldRef() = Zero;

    // Adjust boundary conditions based on the number of steps to be executed
    // and interpolate to faces
    setBoundaryConditions();

    // Solve the elasticity equations in a stepped manner
    for (label istep = 0; istep < nSteps_; ++istep)
    {
        Info<< "Step " << istep << endl;

        // Update diffusivity
        diffusivityPtr_->correct();
        const surfaceScalarField E(diffusivityPtr_->operator()());
        const surfaceVectorField& Sf = fvMesh_.Sf();
        for (label iter = 0; iter < nIters_; ++iter)
        {
            Info<< "Iteration " << iter << endl;
            cellMotionU_.storePrevIter();
            fvVectorMatrix dEqn
            (
                fvm::laplacian(2*E, cellMotionU_)
              + fvc::div(2*E*(fvc::interpolate(fvc::grad(cellMotionU_)) & Sf))
              - fvc::div(E*fvc::interpolate(fvc::div(cellMotionU_))*Sf)
            );

            scalar residual = mag(dEqn.solve().initialResidual());
            cellMotionU_.relax();

            // Print execution time
            fvMesh_.time().printExecutionTime(Info);

            // Check convergence
            if (residual < tolerance_)
            {
                Info<< "\n***Reached mesh movement convergence limit for step "
                    << istep
                    << " iteration " << iter << "***\n\n";
                break;
            }
        }

        interpolationPtr_->interpolate
        (
            cellMotionU_,
            pointMotionU_
        );

        tmp<vectorField> newPoints
        (
            fvMesh_.points() + pointMotionU_.primitiveField()
        );

        /*
        // Interpolate from cells to points
        interpolationPtr_->interpolate(cellMotionU_, pointMotionU_);

        syncTools::syncPointList
        (
            fvMesh_,
            pointMotionU_.primitiveFieldRef(),
            maxEqOp<vector>(),
            vector::zero
        );

        vectorField newPoints
        (
            mesh().points() + pointMotionU_.primitiveFieldRef()
        );
        */

        // Move points and check mesh
        fvMesh_.movePoints(newPoints);
        fvMesh_.checkMesh(true);

        if (debug)
        {
            Info<< "  Writing new mesh points  " << endl;
            pointIOField points
            (
                IOobject
                (
                    "points",
                    mesh().pointsInstance(),
                    polyMesh::meshSubDir,
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh().points()
            );
            points.write();
        }
    }
}


void Foam::elasticityMotionSolver::movePoints(const pointField&)
{
    // Do nothing
}


void Foam::elasticityMotionSolver::updateMesh(const mapPolyMesh&)
{
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
