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

#include "boundaryAdjointContributionIncompressible.H"
#include "adjointRASModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundaryAdjointContributionIncompressible, 0);
addToRunTimeSelectionTable
(
    boundaryAdjointContribution,
    boundaryAdjointContributionIncompressible,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryAdjointContributionIncompressible::
boundaryAdjointContributionIncompressible
(
    const word& managerName,
    const word& adjointSolverName,
    const word& simulationType,
    const fvPatch& patch
)
:
    boundaryAdjointContribution
    (
        managerName,
        adjointSolverName,
        simulationType,
        patch
    ),
    objectiveManager_
    (
        patch_.patch().boundaryMesh().mesh().
            lookupObjectRef<objectiveManager>(managerName)
    ),
    primalVars_
    (
        patch_.patch().boundaryMesh().mesh().
            lookupObject<incompressibleAdjointSolver>(adjointSolverName).
                getPrimalVars()
    ),
    adjointSolver_
    (
        patch_.patch().boundaryMesh().mesh().
            lookupObject<incompressibleAdjointSolver>(adjointSolverName)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> boundaryAdjointContributionIncompressible::velocitySource()
{
    // Objective function contribution
    tmp<vectorField> tsource =
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdv,
            &objectiveIncompressible::hasBoundarydJdv
        );
    vectorField& source = tsource.ref();

    // Turbulence model differentiation contribution.
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();
    source += adjointRAS().adjointMomentumBCSource()[patch_.index()];

    return tsource;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::pressureSource()
{
    // Objective function contribution
    tmp<scalarField> tsource =
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdvn,
            &objectiveIncompressible::hasBoundarydJdvn
        );

    scalarField& source = tsource.ref();

    // Turbulence model differentiation contribution.
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();
    const vectorField& adjointTurbulenceContr =
        adjointRAS().adjointMomentumBCSource()[patch_.index()];

    source += adjointTurbulenceContr & patch_.nf();

    return (tsource);
}


tmp<vectorField>
boundaryAdjointContributionIncompressible::tangentVelocitySource()
{
    // Objective function contribution
    tmp<vectorField> tsource =
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdvt,
            &objectiveIncompressible::hasBoundarydJdvt
        );

    vectorField& source = tsource.ref();

    // Turbulence model differentiation contribution.
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();
    const vectorField& adjointTurbulenceContr =
        adjointRAS().adjointMomentumBCSource()[patch_.index()];

    tmp<vectorField> tnf = patch_.nf();
    const vectorField& nf = tnf();

    source += adjointTurbulenceContr - (adjointTurbulenceContr & nf)*nf;

    return (tsource);
}


tmp<vectorField>
boundaryAdjointContributionIncompressible::normalVelocitySource()
{
    return
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdp,
            &objectiveIncompressible::hasBoundarydJdp
        );

}


tmp<scalarField> boundaryAdjointContributionIncompressible::energySource()
{
    return
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdT,
            &objectiveIncompressible::hasBoundarydJdT
        );

}


tmp<scalarField>
boundaryAdjointContributionIncompressible::adjointTMVariable1Source()
{
    return
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdTMvar1,
            &objectiveIncompressible::hasBoundarydJdTMVar1
        );

}


tmp<scalarField>
boundaryAdjointContributionIncompressible::adjointTMVariable2Source()
{
    return
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdTMvar2,
            &objectiveIncompressible::hasBoundarydJdTMVar2
        );

}


tmp<scalarField>
boundaryAdjointContributionIncompressible::dJdnut()
{
    return
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdnut,
            &objectiveIncompressible::hasBoundarydJdnut
        );
}


tmp<tensorField>
boundaryAdjointContributionIncompressible::dJdGradU()
{
    return
        sumContributions
        (
            objectiveManager_.getObjectiveFunctions(),
            &objectiveIncompressible::boundarydJdGradU,
            &objectiveIncompressible::hasBoundarydJdGradU
        );
}


tmp<scalarField> boundaryAdjointContributionIncompressible::momentumDiffusion()
{

    return adjointVars().adjointTurbulence()().nuEff(patch_.index());
}


tmp<scalarField> boundaryAdjointContributionIncompressible::laminarDiffusivity()
{
    auto tnu = tmp<scalarField>::New(patch_.size(), Zero);
    auto& nu = tnu.ref();

    const autoPtr<incompressible::turbulenceModel>& turbulenceModel =
        primalVars_.turbulence();

    nu = turbulenceModel().nu()().boundaryField()[patch_.index()];

    return tnu;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::thermalDiffusion()
{
    /*
    const polyMesh& mesh = patch_.patch().boundaryMesh().mesh();
    const compressible::turbulenceModel& turbulenceModel =
        mesh.lookupObject<compressible::turbulenceModel>("turbulenceModel");
    tmp<scalarField> talphaEff = turbulenceModel.alphaEff(patch_.index());
    */

    WarningInFunction
        << "no abstract thermalDiffusion is implemented. Returning zero field";

    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContributionIncompressible::wallDistance()
{
    return primalVars_.turbulence()->y()[patch_.index()];
}


tmp<scalarField>
boundaryAdjointContributionIncompressible::TMVariable1Diffusion()
{
    return
        adjointVars().adjointTurbulence()->diffusionCoeffVar1(patch_.index());

}


tmp<scalarField>
boundaryAdjointContributionIncompressible::TMVariable2Diffusion()
{
    return
        adjointVars().adjointTurbulence()->diffusionCoeffVar2(patch_.index());
}


tmp<scalarField> boundaryAdjointContributionIncompressible::TMVariable1()
{
    return
        primalVars_.RASModelVariables()->TMVar1().
            boundaryField()[patch_.index()];
}


tmp<scalarField> boundaryAdjointContributionIncompressible::TMVariable2()
{
    return
        primalVars_.RASModelVariables()->TMVar2().
            boundaryField()[patch_.index()];
}


const fvPatchVectorField& boundaryAdjointContributionIncompressible::Ub() const
{
    return primalVars_.U().boundaryField()[patch_.index()];
}


const fvPatchScalarField& boundaryAdjointContributionIncompressible::pb() const
{
    return primalVars_.p().boundaryField()[patch_.index()];
}


const fvsPatchScalarField&
boundaryAdjointContributionIncompressible::phib() const
{
    return primalVars_.phi().boundaryField()[patch_.index()];
}


tmp<fvPatchScalarField>
boundaryAdjointContributionIncompressible::turbulentDiffusivity() const
{
    return primalVars_.RASModelVariables()().nutPatchField(patch_.index());
}


const fvPatchVectorField& boundaryAdjointContributionIncompressible::Uab() const
{
    return adjointVars().UaInst().boundaryField()[patch_.index()];
}


const fvPatchScalarField& boundaryAdjointContributionIncompressible::pab() const
{
    return adjointVars().paInst().boundaryField()[patch_.index()];
}


const fvsPatchScalarField&
boundaryAdjointContributionIncompressible::phiab() const
{
    return adjointVars().phiaInst().boundaryField()[patch_.index()];
}


const word boundaryAdjointContributionIncompressible::primalSolverName() const
{
    return primalVars_.solverName();
}


const word boundaryAdjointContributionIncompressible::adjointSolverName() const
{
    return adjointVars().solverName();
}


const incompressibleVars&
boundaryAdjointContributionIncompressible::primalVars() const
{
    return primalVars_;
}


const incompressibleAdjointVars&
boundaryAdjointContributionIncompressible::adjointVars() const
{
    return adjointSolver_.getAdjointVars();
}


objectiveManager&
boundaryAdjointContributionIncompressible::getObjectiveManager()
{
    return objectiveManager_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
