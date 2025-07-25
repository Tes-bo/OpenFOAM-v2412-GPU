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

#include "adjointMeshMovementSolver.H"
#include "adjointEikonalSolver.H"
#include "adjointSolver.H"
#include "fvc.H"
#include "fvm.H"
#include "ShapeSensitivitiesBase.H"
#include "reverseLinear.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(adjointMeshMovementSolver, 0);

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void adjointMeshMovementSolver::read()
{
    iters_ = dict_.getOrDefault<label>("iters", 1000);
    tolerance_ = dict_.getOrDefault<scalar>("tolerance", 1.e-06);
}


void adjointMeshMovementSolver::setSource()
{
    volTensorField& gradDxDbMult = adjointSensitivity_.gradDxDbMult()();

    // Add part related to the adjoint eikaonal equation, if necessary
    const autoPtr<adjointEikonalSolver>& eikonalSolver =
        adjointSensitivity_.getAdjointEikonalSolver();
    if (eikonalSolver)
    {
        gradDxDbMult += eikonalSolver->getFISensitivityTerm();
    }

    source_ -=
        fvc::div
        (
            mesh_.Sf()
          & reverseLinear<tensor>(mesh_).interpolate(gradDxDbMult)
        );

    // Terms from objectives defined in (part of the) internal field
    PtrList<objective>& functions =
        adjointSensitivity_.getAdjointSolver().getObjectiveManager().
            getObjectiveFunctions();
    for (objective& func : functions)
    {
        if (func.hasDivDxDbMult())
        {
            source_ -= func.weight()*fvc::grad(func.divDxDbMultiplier());
        }
    }

    // Terms from fvOptions
    source_.primitiveFieldRef() += adjointSensitivity_.optionsDxDbMult()();
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

adjointMeshMovementSolver::adjointMeshMovementSolver
(
    const fvMesh& mesh,
    const dictionary& dict,
    ShapeSensitivitiesBase& adjointSensitivity
)
:
    mesh_(mesh),
    dict_(dict.subOrEmptyDict("adjointMeshMovementSolver")),
    meshMovementSensPtr_(createZeroBoundaryPtr<vector>(mesh)),
    adjointSensitivity_(adjointSensitivity),
    ma_
    (
        variablesSet::autoCreateMeshMovementField
        (
            mesh_,
            adjointSensitivity.getAdjointSolver().useSolverNameForFields()
          ? ("ma" + adjointSensitivity.getAdjointSolver().solverName())
          : "ma",
            adjointSensitivity.getAdjointSolver().maDimensions()
        )
    ),
    source_
    (
        IOobject
        (
            "sourceadjointMeshMovement",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            adjointSensitivity.getAdjointSolver().maDimensions()/sqr(dimLength),
            Zero
        )
    ),
    iters_(0),
    tolerance_(Zero)
{
    read();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool adjointMeshMovementSolver::readDict
(
    const dictionary& dict
)
{
    dict_ = dict.subOrEmptyDict("adjointMeshMovementSolver");
    read();

    return true;
}


void adjointMeshMovementSolver::solve()
{
    setSource();

    // Iterate the adjoint to the mesh movement equation
    for (label iter = 0; iter < iters_; iter++)
    {
        Info<< "adjoint Mesh Movement Iteration: " << iter << endl;

        fvVectorMatrix maEqn
        (
            fvm::laplacian(ma_) + source_
        );

        maEqn.boundaryManipulate(ma_.boundaryFieldRef());

        scalar residual =
            mag(Foam::solve(maEqn, mesh_.solverDict("ma")).initialResidual());

        Info<< "Max ma " << gMax(mag(ma_)()) << endl;

        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance_)
        {
            Info<< "\n***Reached adjoint mesh movement convergence limit, "
                   "iteration " << iter << "***\n\n";
            break;
        }
    }
    ma_.write();
}


void adjointMeshMovementSolver::reset()
{
    source_ == dimensionedVector(source_.dimensions(), Zero);
    meshMovementSensPtr_() = Zero;
}


boundaryVectorField& adjointMeshMovementSolver::meshMovementSensitivities()
{
    boundaryVectorField& meshMovementSens = meshMovementSensPtr_();

    for
    (
        const label patchi
      : adjointSensitivity_.geometryVariationIntegrationPatches()
    )
    {
        // No surface area included.
        // Will be added during the assembly of the sensitivities
        meshMovementSens[patchi] = -ma_.boundaryField()[patchi].snGrad();
    }

    return meshMovementSens;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
