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

#include "laplacianMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laplacianMotionSolver, 1);

    addToRunTimeSelectionTable
    (
        motionSolver,
        laplacianMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laplacianMotionSolver::laplacianMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
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
    nIters_(this->coeffDict().get<label>("iters")),
    tolerance_(this->coeffDict().get<scalar>("tolerance"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::laplacianMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellMotionU_,
        pointMotionU_
    );

    syncTools::syncPointList
    (
        fvMesh_,
        pointMotionU_.primitiveFieldRef(),
        maxEqOp<vector>(),
        vector::zero
    );

    tmp<vectorField> tcurPoints
    (
        fvMesh_.points() + pointMotionU_.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::laplacianMotionSolver::solve()
{
    diffusivityPtr_->correct();

    // Iteratively solve the Laplace equation, to account for non-orthogonality
    for (label iter = 0; iter < nIters_; ++iter)
    {
        Info<< "Iteration " << iter << endl;
        fvVectorMatrix dEqn
        (
            fvm::laplacian
            (
                dimensionedScalar("viscosity", dimViscosity, 1.0)
              * diffusivityPtr_->operator()(),
                cellMotionU_,
                "laplacian(diffusivity,cellMotionU)"
            )
        );

        scalar residual = mag(dEqn.solve().initialResidual());

        // Print execution time
        fvMesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance_)
        {
            Info<< "\n***Reached mesh movement convergence limit at"
                << " iteration " << iter << "***\n\n";
            break;
        }
    }
}


void Foam::laplacianMotionSolver::setBoundaryConditions()
{
    pointMotionU_.boundaryFieldRef().updateCoeffs();
    auto& cellMotionUbf = cellMotionU_.boundaryFieldRef();

    forAll(cellMotionU_.boundaryField(), pI)
    {
        fvPatchVectorField& bField = cellMotionUbf[pI];
        if (isA<fixedValueFvPatchVectorField>(bField))
        {
            const pointField& points = fvMesh_.points();
            const polyPatch& patch = fvMesh_.boundaryMesh()[pI];
            forAll(bField, fI)
            {
                bField[fI] = patch[fI].average(points, pointMotionU_);
            }
        }
    }
}


void Foam::laplacianMotionSolver::movePoints(const pointField&)
{
    // Do nothing
}


void Foam::laplacianMotionSolver::updateMesh(const mapPolyMesh&)
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
