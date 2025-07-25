/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "objectiveIncompressible.H"
#include "incompressiblePrimalSolver.H"
#include "incompressibleAdjointSolver.H"
#include "createZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveIncompressible, 0);
defineRunTimeSelectionTable(objectiveIncompressible, dictionary);
addToRunTimeSelectionTable
(
    objective,
    objectiveIncompressible,
    objective
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveIncompressible::objectiveIncompressible
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objective(mesh, dict, adjointSolverName, primalSolverName),

    vars_
    (
        mesh.lookupObject<incompressiblePrimalSolver>(primalSolverName).
            getIncoVars()
    ),

    // Initialize pointers to nullptr.
    // Not all of them are required for each objective function.
    // Each child should allocate whatever is needed.

    // Field adjoint Eqs
    dJdvPtr_(nullptr),
    dJdpPtr_(nullptr),
    dJdTPtr_(nullptr),
    dJdTMvar1Ptr_(nullptr),
    dJdTMvar2Ptr_(nullptr),

    // Adjoint boundary conditions
    bdJdvPtr_(nullptr),
    bdJdvnPtr_(nullptr),
    bdJdvtPtr_(nullptr),
    bdJdpPtr_(nullptr),
    bdJdTPtr_(nullptr),
    bdJdTMvar1Ptr_(nullptr),
    bdJdTMvar2Ptr_(nullptr),
    bdJdnutPtr_(nullptr),
    bdJdGradUPtr_(nullptr)
{
    computeMeanFields_ = vars_.computeMeanFields();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<objectiveIncompressible> objectiveIncompressible::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Creating objective function : " << dict.dictName()
        << " of type " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "objectiveIncompressible",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<objectiveIncompressible>
    (
        ctorPtr(mesh, dict, adjointSolverName, primalSolverName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void objectiveIncompressible::doNormalization()
{
    if (normalize_ && normFactor_)
    {
        const scalar oneOverNorm(1./normFactor_());

        if (hasdJdv())
        {
            dJdvPtr_().primitiveFieldRef() *= oneOverNorm;
        }
        if (hasdJdp())
        {
            dJdpPtr_().primitiveFieldRef() *= oneOverNorm;
        }
        if (hasdJdT())
        {
            dJdTPtr_().primitiveFieldRef() *= oneOverNorm;
        }
        if (hasdJdTMVar1())
        {
            dJdTMvar1Ptr_().primitiveFieldRef() *= oneOverNorm;
        }
        if (hasdJdTMVar2())
        {
            dJdTMvar2Ptr_().primitiveFieldRef() *= oneOverNorm;
        }
        if (hasBoundarydJdv())
        {
            bdJdvPtr_() *= oneOverNorm;
        }
        if (hasBoundarydJdvn())
        {
            bdJdvnPtr_() *= oneOverNorm;
        }
        if (hasBoundarydJdvt())
        {
            bdJdvtPtr_() *= oneOverNorm;
        }
        if (hasBoundarydJdp())
        {
            bdJdpPtr_() *= oneOverNorm;
        }
        if (hasBoundarydJdT())
        {
            bdJdTPtr_() *= oneOverNorm;
        }
        if (hasBoundarydJdTMVar1())
        {
            bdJdTMvar1Ptr_() *= oneOverNorm;
        }
        if (hasBoundarydJdTMVar2())
        {
            bdJdTMvar2Ptr_() *= oneOverNorm;
        }
        if (hasBoundarydJdnut())
        {
            bdJdnutPtr_() *= oneOverNorm;
        }
        if (hasBoundarydJdGradU())
        {
            bdJdGradUPtr_() *= oneOverNorm;
        }

        // Normalize geometric fields
        objective::doNormalization();
    }
}


void objectiveIncompressible::update()
{
    // Update geometric fields
    objective::update();

    // Update mean values here since they might be used in the
    // subsequent functions
    update_meanValues();

    // volFields
    update_dJdv();
    update_dJdp();
    update_dJdT();
    update_dJdTMvar1();
    update_dJdTMvar2();

    // boundaryFields
    update_boundarydJdv();
    update_boundarydJdvn();
    update_boundarydJdvt();
    update_boundarydJdp();
    update_boundarydJdT();
    update_boundarydJdTMvar1();
    update_boundarydJdTMvar2();
    update_boundarydJdnut();
    update_boundarydJdGradU();

    // Divide everything with normalization factor
    doNormalization();

    // Set objective as not computed, for the next optimisation cycle
    computed_ = false;
}


void objectiveIncompressible::nullify()
{
    if (!nullified_)
    {
        if (hasdJdv())
        {
            dJdvPtr_() == dimensionedVector(dJdvPtr_().dimensions(), Zero);
        }
        if (hasdJdp())
        {
            dJdpPtr_() == dimensionedScalar(dJdpPtr_().dimensions(), Zero);
        }
        if (hasdJdT())
        {
            dJdTPtr_() == dimensionedScalar(dJdTPtr_().dimensions(), Zero);
        }
        if (hasdJdTMVar1())
        {
            dJdTMvar1Ptr_() ==
                dimensionedScalar(dJdTMvar1Ptr_().dimensions(), Zero);
        }
        if (hasdJdTMVar2())
        {
            dJdTMvar2Ptr_() ==
                dimensionedScalar(dJdTMvar2Ptr_().dimensions(), Zero);
        }
        if (hasBoundarydJdv())
        {
            bdJdvPtr_() == vector::zero;
        }
        if (hasBoundarydJdvn())
        {
            bdJdvnPtr_() == scalar(0);
        }
        if (hasBoundarydJdvt())
        {
            bdJdvtPtr_() == vector::zero;
        }
        if (hasBoundarydJdp())
        {
            bdJdpPtr_() == vector::zero;
        }
        if (hasBoundarydJdT())
        {
            bdJdTPtr_() == scalar(0);
        }
        if (hasBoundarydJdTMVar1())
        {
            bdJdTMvar1Ptr_() == scalar(0);
        }
        if (hasBoundarydJdTMVar2())
        {
            bdJdTMvar2Ptr_() == scalar(0);
        }
        if (hasBoundarydJdnut())
        {
            bdJdnutPtr_() == scalar(0);
        }
        if (hasBoundarydJdGradU())
        {
            bdJdGradUPtr_() == tensor::zero;
        }

        // Nullify geometric fields and sets nullified_ to true
        objective::nullify();
    }
}


void objectiveIncompressible::allocatedJdTurbulence()
{
    // Figure out the availability of the adjoint turbulence model variables
    // through their primal counterparts, since the contructor of the adjoint
    // solver has not been completed yet
    const incompressiblePrimalSolver& primSolver =
        mesh_.lookupObject<incompressiblePrimalSolver>(primalSolverName_);
    const autoPtr<incompressible::RASModelVariables>& rasVars =
        primSolver.getIncoVars().RASModelVariables();

    if (rasVars().hasTMVar1())
    {
        const dimensionSet primalVarDims = rasVars->TMVar1Inst().dimensions();
        const dimensionSet sourceDims = dimArea/pow3(dimTime)/primalVarDims;
        dJdTMvar1Ptr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                "ATMSource1",
                sourceDims
            )
        );
    }
    if (rasVars().hasTMVar2())
    {
        const dimensionSet primalVarDims = rasVars->TMVar2Inst().dimensions();
        const dimensionSet sourceDims = dimArea/pow3(dimTime)/primalVarDims;
        dJdTMvar2Ptr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                "ATMSource2",
                sourceDims
            )
        );
    }
}


void objectiveIncompressible::checkCellZonesSize
(
    const labelList& zoneIDs
) const
{
    label nCells(0);
    for (const label zI : zoneIDs)
    {
        nCells += mesh_.cellZones()[zI].size();
    }
    reduce(nCells, sumOp<label>());
    if (!nCells)
    {
        FatalErrorInFunction
            << "Provided cellZones include no cells"
            << exit(FatalError);
    }
}


void objectiveIncompressible::update_dJdTMvar
(
    autoPtr<volScalarField>& dJdTMvarPtr,
    tmp<volScalarField>
        (incompressibleAdjoint::adjointRASModel::*JacobianFunc)() const,
    const volScalarField& JacobianMultiplier,
    const labelList& zones
)
{
    if (dJdTMvarPtr.good())
    {
        // nut Jacobians are currently computed in the adjoint turbulence
        // models, though they would be better placed within the primal
        // turbulence model.
        // Safeguard the computation until the adjoint solver is complete
        if (mesh_.foundObject<incompressibleAdjointSolver>(adjointSolverName_))
        {
            const incompressibleAdjointSolver& adjSolver =
                mesh_.lookupObject<incompressibleAdjointSolver>
                    (adjointSolverName_);
            const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
                adjSolver.getAdjointVars().adjointTurbulence();

            tmp<volScalarField> tnutJacobian = (adjointRAS->*JacobianFunc)();
            const volScalarField& nutJacobian = tnutJacobian();

            volScalarField& dJdTMvar = dJdTMvarPtr();

            for (const label zI : zones)
            {
                const cellZone& zoneI = mesh_.cellZones()[zI];
                for (const label cellI : zoneI)
                {
                    dJdTMvar[cellI] =
                        JacobianMultiplier[cellI]*nutJacobian[cellI];
                }
            }
        }
        else
        {
            WarningInFunction
                << "Skipping the computation of nutJacobian until "
                << "the adjoint solver is complete"
                << endl;
        }
    }
}


void objectiveIncompressible::addSource(fvVectorMatrix& matrix)
{
    if (fieldNames_.found(matrix.psi().name()) && hasdJdv())
    {
        matrix += weight()*dJdv();
    }
}


bool objectiveIncompressible::write(const bool valid) const
{
    return objective::write(valid);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
