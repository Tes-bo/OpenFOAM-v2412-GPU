/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "objectivePowerDissipation.H"
#include "incompressibleAdjointSolver.H"
#include "createZeroField.H"
#include "topOVariablesBase.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectivePowerDissipation, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectivePowerDissipation,
    dictionary
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void objectivePowerDissipation::populateFieldNames()
{
    if (fieldNames_.size() == 1)
    {
        const incompressibleAdjointSolver& adjSolver =
            mesh_.lookupObject<incompressibleAdjointSolver>(adjointSolverName_);
        const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
            adjSolver.getAdjointVars().adjointTurbulence();
        const wordList& baseNames =
            adjointRAS().getAdjointTMVariablesBaseNames();
        forAll(baseNames, nI)
        {
            fieldNames_.push_back
                (adjSolver.extendedVariableName(baseNames[nI]));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectivePowerDissipation::objectivePowerDissipation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    zones_(mesh_.cellZones().indices(dict.get<wordRes>("zones")))
{
    // Append Ua name to fieldNames
    fieldNames_.setSize
    (
        1,
        mesh_.lookupObject<solver>(adjointSolverName_).
            extendedVariableName("Ua")
    );

    // Check if cellZones provided include at least one cell
    checkCellZonesSize(zones_);

    // Allocate dJdTMvar1Ptr_ and dJdTMvar2Ptr_ if needed
    allocatedJdTurbulence();

    // Allocate source term to the adjoint momentum equations
    dJdvPtr_.reset
    (
        createZeroFieldPtr<vector>
        (
            mesh_,
            "dJdv" + type(),
            dimLength/sqr(dimTime)
        )
    );
    // Allocate terms to be added to volume-based sensitivity derivatives
    divDxDbMultPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "divDxDbMult" + objectiveName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    gradDxDbMultPtr_.reset
    (
        createZeroFieldPtr<tensor>
        (
            mesh_,
            ("gradDxdbMult" + type()),
            sqr(dimLength)/pow3(dimTime)
        )
    );

    // Allocate direct sensitivity contributions for topology optimisation
    dJdbPtr_.reset(createZeroFieldPtr<scalar>(mesh_, "dJdb", dimless));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectivePowerDissipation::J()
{
    J_ = Zero;

    // References
    const volVectorField& U = vars_.UInst();
    const autoPtr<incompressible::turbulenceModel>& turb = vars_.turbulence();
    const scalarField& V = mesh_.V().field();

    volScalarField integrand(turb->nuEff()*magSqr(twoSymm(fvc::grad(U))));

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        scalarField VZone(V, zoneI);
        scalarField integrandZone(integrand.primitiveField(), zoneI);

        J_ += 0.5*gSum(integrandZone*VZone);
        if (mesh_.foundObject<topOVariablesBase>("topOVars"))
        {
            const topOVariablesBase& vars =
                mesh_.lookupObject<topOVariablesBase>("topOVars");
            const volScalarField& beta = vars.beta();
            scalar porosityContr = Zero;
            for (const label cellI : zoneI)
            {
                porosityContr += beta[cellI]*magSqr(U[cellI])*V[cellI];
            }
            porosityContr *= vars.getBetaMax();
            J_ += returnReduce(porosityContr, sumOp<scalar>());
        }
    }

    return J_;
}


void objectivePowerDissipation::update_dJdv()
{
    dJdvPtr_().primitiveFieldRef() = Zero;

    const volVectorField& U = vars_.U();
    const autoPtr<incompressible::turbulenceModel>& turb = vars_.turbulence();
    tmp<volSymmTensorField> S = twoSymm(fvc::grad(U));

    // Add source from possible dependencies of nut on U
    if (mesh_.foundObject<incompressibleAdjointSolver>(adjointSolverName_))
    {
        const incompressibleAdjointSolver& adjSolver =
            mesh_.lookupObject<incompressibleAdjointSolver>
                (adjointSolverName_);
        const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
            adjSolver.getAdjointVars().adjointTurbulence();
        tmp<volScalarField> dnutdUMult = 0.5*magSqr(S());
        tmp<volVectorField> dnutdU = adjointRAS->nutJacobianU(dnutdUMult);
        if (dnutdU)
        {
            for (const label zI : zones_)
            {
                const cellZone& zoneI = mesh_.cellZones()[zI];
                for (const label cellI : zoneI)
                {
                    dJdvPtr_()[cellI] = dnutdU()[cellI];
                }
            }
        }
    }

    // Add source from the strain rate magnitude
    volVectorField integrand(-2.0*fvc::div(turb->nuEff()*S));
    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            dJdvPtr_()[cellI] += integrand[cellI];
        }
    }

    // Add source from porosity dependencies
    if (mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        const topOVariablesBase& vars =
            mesh_.lookupObject<topOVariablesBase>("topOVars");
        const volScalarField& beta = vars.beta();
        const scalar betaMax = vars.getBetaMax();
        for (const label zI : zones_)
        {
            const cellZone& zoneI = mesh_.cellZones()[zI];
            for (const label cellI : zoneI)
            {
                dJdvPtr_()[cellI] += 2*betaMax*beta[cellI]*U[cellI];
            }
        }
    }
}


void objectivePowerDissipation::update_dJdTMvar1()
{
    const volVectorField& U = vars_.U();
    volScalarField JacobianMultiplier(0.5*magSqr(twoSymm(fvc::grad(U))));
    update_dJdTMvar
    (
        dJdTMvar1Ptr_,
        &incompressibleAdjoint::adjointRASModel::nutJacobianTMVar1,
        JacobianMultiplier,
        zones_
    );
}


void objectivePowerDissipation::update_dJdTMvar2()
{
    const volVectorField& U = vars_.U();
    volScalarField JacobianMultiplier(0.5*magSqr(twoSymm(fvc::grad(U))));
    update_dJdTMvar
    (
        dJdTMvar2Ptr_,
        &incompressibleAdjoint::adjointRASModel::nutJacobianTMVar2,
        JacobianMultiplier,
        zones_
    );
}


void objectivePowerDissipation::update_divDxDbMultiplier()
{
    // References
    volScalarField& divDxDbMult = divDxDbMultPtr_();
    const volVectorField& U = vars_.U();
    const autoPtr<incompressible::turbulenceModel>& turb = vars_.turbulence();
    volScalarField integrand(0.5*turb->nuEff()*magSqr(twoSymm(fvc::grad(U))));

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            divDxDbMult[cellI] = integrand[cellI];
        }
    }
    divDxDbMult.correctBoundaryConditions();
}


void objectivePowerDissipation::update_gradDxDbMultiplier()
{
    // References
    volTensorField& gradDxDbMult = gradDxDbMultPtr_();
    const volVectorField& U = vars_.U();
    const autoPtr<incompressible::turbulenceModel>& turb = vars_.turbulence();
    tmp<volTensorField> gradU = fvc::grad(U);
    volTensorField integrand(-2.0*turb->nuEff()*(gradU() & twoSymm(gradU())));

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            gradDxDbMult[cellI] = integrand[cellI];
        }
    }
    gradDxDbMult.correctBoundaryConditions();
}


void objectivePowerDissipation::update_dJdb()
{
    if (mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        scalarField& dJdb = dJdbPtr_().primitiveFieldRef();
        dJdb = Zero;
        const volVectorField& U = vars_.UInst();
        const topOVariablesBase& vars =
            mesh_.lookupObject<topOVariablesBase>("topOVars");
        const scalar betaMax = vars.getBetaMax();
        for (const label zI : zones_)
        {
            const cellZone& zoneI = mesh_.cellZones()[zI];
            for (const label cellI : zoneI)
            {
                dJdb[cellI] += betaMax*magSqr(U[cellI]);
            }
        }
    }
}




void objectivePowerDissipation::addSource(fvScalarMatrix& matrix)
{
    populateFieldNames();
    const label fieldI = fieldNames_.find(matrix.psi().name());
    if (fieldI == 1)
    {
        matrix += weight()*dJdTMvar1Ptr_();
    }
    if (fieldI == 2)
    {
        matrix += weight()*dJdTMvar2Ptr_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
