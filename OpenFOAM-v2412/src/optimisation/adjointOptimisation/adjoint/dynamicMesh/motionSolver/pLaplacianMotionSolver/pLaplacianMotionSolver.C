/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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

#include "pLaplacianMotionSolver.H"
#include "motionInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pLaplacianMotionSolver, 1);

    addToRunTimeSelectionTable
    (
        motionSolver,
        pLaplacianMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pLaplacianMotionSolver::pLaplacianMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    useFixedValuePointMotionUBCs_
        (coeffDict().getOrDefault<bool>("useFixedValuePointMotionUBCs", false)),
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
        word
        (
            useFixedValuePointMotionUBCs_
          ? fixedValuePointPatchVectorField::typeName
          : calculatedPointPatchField<vector>::typeName
        )
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
    nIters_(this->coeffDict().get<label>("iters")),
    tolerance_(this->coeffDict().get<scalar>("tolerance")),
    toleranceIntermediate_
    (
        this->coeffDict().
            getOrDefault<scalar>("toleranceIntermediate", 100*tolerance_)
    ),
    exponent_(this->coeffDict().get<label>("exponent"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pLaplacianMotionSolver::setBoundaryConditions()
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



Foam::tmp<Foam::pointField> Foam::pLaplacianMotionSolver::curPoints() const
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
        fvMesh_.points() + pointMotionU_.internalField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::pLaplacianMotionSolver::solve()
{
//  setBoundaryConditions();

    for (label exp = 2; exp < exponent_ + 1; ++exp)
    {
        scalar tolerance
            (exp == exponent_ ? tolerance_ : toleranceIntermediate_);
        Info<< "Solving for exponent " << exp << endl;

        for (label iter = 0; iter < nIters_; ++iter)
        {
            Info<< "Iteration " << iter << endl;
            cellMotionU_.storePrevIter();
            volScalarField gamma(pow(mag(fvc::grad(cellMotionU_)), exp - 2));
            gamma.correctBoundaryConditions();
            fvVectorMatrix dEqn
            (
                fvm::laplacian(gamma, cellMotionU_)
            );

            scalar residual = mag(dEqn.solve().initialResidual());

            cellMotionU_.relax();

            // Print execution time
            fvMesh_.time().printExecutionTime(Info);

            // Check convergence
            if (residual < tolerance)
            {
                Info<< "\n***Reached mesh movement convergence limit at"
                    << " iteration " << iter << "***\n\n";
                if (debug)
                {
                    gamma.write();
                }
                break;
            }
        }
    }
}


void Foam::pLaplacianMotionSolver::movePoints(const pointField&)
{
    // Do nothing
}


void Foam::pLaplacianMotionSolver::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


// ************************************************************************* //
