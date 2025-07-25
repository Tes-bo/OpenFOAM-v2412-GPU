/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 PCOpt/NTUA
    Copyright (C) 2023 FOSS GP
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

#include "adjointNull.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointNull, 0);
    addToRunTimeSelectionTable
    (
        adjointSolver,
        adjointNull,
        adjointSolver
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::adjointNull::preCalculateSensitivities()
{
    adjointSensitivity_->accumulateIntegrand(scalar(1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointNull::adjointNull
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
:
    adjointSolver
    (
        mesh,
        managerType,
        dict,
        primalSolverName,
        solverName
    )
{
    allocateSensitivities();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::adjointNull> Foam::adjointNull::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
{
    return autoPtr<adjointNull>
    (
        new adjointNull
        (
            mesh,
            managerType,
            dict,
            primalSolverName,
            solverName
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::adjointNull::simulationType() const
{
    return word("null");
}


Foam::dimensionSet Foam::adjointNull::maDimensions() const
{
    return pow3(dimLength/dimTime);
}


void Foam::adjointNull::solveIter()
{
    // Does nothing
}


void Foam::adjointNull::solve()
{
    // Does nothing
}


bool Foam::adjointNull::loop()
{
    return false;
}


void Foam::adjointNull::updatePrimalBasedQuantities()
{
    // Update objective function related quantities
    objectiveManager_.updateAndWrite();
}


void Foam::adjointNull::accumulateGradDxDbMultiplier
(
    volTensorField& gradDxDbMult,
    const scalar dt
)
{
    tmp<volTensorField> tsens
    (
        tmp<volTensorField>::New
        (
            IOobject
            (
                "flowTerm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(sqr(dimLength)/pow3(dimTime), Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    volTensorField& sens = tsens.ref();

    PtrList<objective>& functions = objectiveManager_.getObjectiveFunctions();

    for (objective& objI : functions)
    {
        if (objI.hasGradDxDbMult())
        {
            sens += objI.weight()*objI.gradDxDbMultiplier();
        }
    }
    sens.correctBoundaryConditions();

    gradDxDbMult += sens.T()*dt;
}


void Foam::adjointNull::accumulateDivDxDbMultiplier
(
    autoPtr<scalarField>& divDxDbMult,
    const scalar dt
)
{
    PtrList<objective>& functions = objectiveManager_.getObjectiveFunctions();
    for (objective& func : functions)
    {
        if (func.hasDivDxDbMult())
        {
            divDxDbMult() +=
                func.weight()*func.divDxDbMultiplier().primitiveField()*dt;
        }
    }
}


void Foam::adjointNull::accumulateGeometryVariationsMultipliers
(
    autoPtr<boundaryVectorField>& dSfdbMult,
    autoPtr<boundaryVectorField>& dnfdbMult,
    autoPtr<boundaryVectorField>& dxdbDirectMult,
    autoPtr<pointBoundaryVectorField>& pointDxDbDirectMult,
    const labelHashSet& sensitivityPatchIDs,
    const scalar dt
)
{
    PtrList<objective>& functions = objectiveManager_.getObjectiveFunctions();
    for (const label patchI : sensitivityPatchIDs)
    {
        const scalarField magSfDt(mesh_.boundary()[patchI].magSf()*dt);
        for (objective& func : functions)
        {
            const scalar wei(func.weight());
            if (func.hasdSdbMult())
            {
                dSfdbMult()[patchI] += wei*func.dSdbMultiplier(patchI)*dt;
            }
            if (func.hasdndbMult())
            {
                dnfdbMult()[patchI] += wei*func.dndbMultiplier(patchI)*magSfDt;
            }
            if (func.hasdxdbDirectMult())
            {
                dxdbDirectMult()[patchI] +=
                    wei*func.dxdbDirectMultiplier(patchI)*magSfDt;
            }
        }
    }
}


void Foam::adjointNull::topOSensMultiplier
(
    scalarField& betaMult,
    const word& designVariablesName,
    const scalar dt
)
{
    // Terms resulting directly from the objective function
    PtrList<objective>& functions = objectiveManager_.getObjectiveFunctions();
    for (objective& objI : functions)
    {
        const scalar weight = objI.weight();
        if (objI.hasdJdb())
        {
            betaMult += weight*objI.dJdb()*dt;
        }

        if (objI.hasdJdbField())
        {
            SubField<scalar> betaSens(objI.dJdbField(), mesh_.nCells(), 0);
            betaMult += weight*betaSens*dt;
        }
    }
}


// ************************************************************************* //
