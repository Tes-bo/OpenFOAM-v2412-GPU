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

#include "incompressibleAdjointSolver.H"
#include "incompressiblePrimalSolver.H"
#include "wallFvPatch.H"
#include "sensitivityTopO.H"
#include "adjointBoundaryConditions.H"
#include "adjointBoundaryConditionsFwd.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleAdjointSolver, 0);
    defineRunTimeSelectionTable(incompressibleAdjointSolver, dictionary);
    addToRunTimeSelectionTable
    (
        adjointSolver,
        incompressibleAdjointSolver,
        adjointSolver
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::incompressibleAdjointSolver::hasBCdxdbMult
(
    const labelHashSet& sensitivityPatchIDs
)
{
    if (hasBCdxdbMult_.bad())
    {
        const volVectorField& Ua = getAdjointVars().Ua();
        hasBCdxdbMult_ = false;
        for (const label patchI : sensitivityPatchIDs)
        {
            const fvPatchVectorField& Uab = Ua.boundaryField()[patchI];
            if (isA<adjointVectorBoundaryCondition>(Uab))
            {
                adjointVectorBoundaryCondition& abc =
                    refCast<adjointVectorBoundaryCondition>
                        (const_cast<fvPatchVectorField&>(Uab));
                tmp<tensorField> dxdbMult = abc.dxdbMult();
                if (dxdbMult)
                {
                    hasBCdxdbMult_ = true;
                    break;
                }
            }
        }
    }
    return hasBCdxdbMult_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleAdjointSolver::incompressibleAdjointSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
:
    adjointSolver(mesh, managerType, dict, primalSolverName, solverName),
    primalVars_
    (
        mesh.lookupObjectRef<incompressiblePrimalSolver>(primalSolverName).
            getIncoVars()
    ),
    ATCModel_(nullptr),
    hasBCdxdbMult_(Switch::INVALID)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::incompressibleAdjointSolver>
Foam::incompressibleAdjointSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
{
    const word solverType(dict.get<word>("solver"));
    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "incompressibleAdjointSolver",
            solverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return
        autoPtr<incompressibleAdjointSolver>
        (
            ctorPtr(mesh, managerType, dict, primalSolverName, solverName)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::incompressibleVars&
Foam::incompressibleAdjointSolver::getPrimalVars() const
{
    return primalVars_;
}



const Foam::incompressibleAdjointVars&
Foam::incompressibleAdjointSolver::getAdjointVars() const
{
    const incompressibleAdjointVars& adjointVars =
        refCast<incompressibleAdjointVars>(const_cast<variablesSet&>(vars_()));
    return adjointVars;
}


Foam::incompressibleAdjointVars&
Foam::incompressibleAdjointSolver::getAdjointVars()
{
    incompressibleAdjointVars& adjointVars =
        refCast<incompressibleAdjointVars>(const_cast<variablesSet&>(vars_()));
    return adjointVars;
}



const Foam::autoPtr<Foam::ATCModel>&
Foam::incompressibleAdjointSolver::getATCModel() const
{
    return ATCModel_;
}


bool Foam::incompressibleAdjointSolver::includeDistance() const
{
    return getAdjointVars().adjointTurbulence()->includeDistance();
}


Foam::dimensionSet Foam::incompressibleAdjointSolver::daDimensions() const
{
    return sqr(dimLength)/pow3(dimTime);
}


Foam::dimensionSet Foam::incompressibleAdjointSolver::maDimensions() const
{
    return pow3(dimLength/dimTime);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleAdjointSolver::adjointEikonalSource()
{
    return getAdjointVars().adjointTurbulence()->distanceSensitivities();
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleAdjointSolver::yWall() const
{
    return getPrimalVars().RASModelVariables()->d();
}


Foam::autoPtr<Foam::ATCModel>& Foam::incompressibleAdjointSolver::getATCModel()
{
    return ATCModel_;
}


void Foam::incompressibleAdjointSolver::updatePrimalBasedQuantities()
{
    if (vars_)
    {
        getAdjointVars().adjointTurbulence()->setChangedPrimalSolution();
        ATCModel_().updatePrimalBasedQuantities();
        getAdjointVars().updatePrimalBasedQuantities();
    }
}


void Foam::incompressibleAdjointSolver::accumulateGradDxDbMultiplier
(
    volTensorField& gradDxDbMult,
    const scalar dt
)
{
    /*
    addProfiling
    (
        incompressibleAdjointSolver,
        "incompressibleAdjointSolver::accumulateGradDxDbMultiplier"
    );
    */
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS
    (
        getAdjointVars().adjointTurbulence()
    );

    const volScalarField& p = primalVars_.p();
    const volVectorField& U = primalVars_.U();
    const volScalarField& pa = getAdjointVars().pa();
    const volVectorField& Ua = getAdjointVars().Ua();

    // We only need to modify the boundaryField of gradU locally.
    // If grad(U) is cached then
    // a. The .ref() call fails since the tmp is initialised from a
    //    const ref
    // b. we would be changing grad(U) for all other places in the code
    //    that need it
    // So, always allocate new memory and avoid registering the new field
    tmp<volTensorField> tgradU =
        volTensorField::New("gradULocal", fvc::grad(U));
    volTensorField& gradU = tgradU.ref();
    volTensorField::Boundary& gradUbf = gradU.boundaryFieldRef();

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
            gradUbf[patchI] = tnf*U.boundaryField()[patchI].snGrad();
        }
    }

    tmp<volScalarField> tnuEff = adjointRAS->nuEff();
    tmp<volSymmTensorField> stress = tnuEff()*twoSymm(gradU);
    tmp<volTensorField> tgradUa = fvc::grad(Ua);

    // Return field
    auto tflowTerm
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
    volTensorField& flowTerm = tflowTerm.ref();
    // Term 3
    flowTerm = - tnuEff*(gradU & twoSymm(tgradUa()));
    // Term 4
    flowTerm += fvc::grad(Ua & stress()) - (tgradUa & stress());

    // Release memory
    stress.clear();

    // Compute dxdb multiplier
    flowTerm +=
        // Term 1, ATC
        ATCModel_->getFISensitivityTerm()
        // Term 2
      - fvc::grad(p)*Ua;

    // Term 5
    flowTerm += pa*tgradU;

    // Term 6, from the adjoint turbulence model
    flowTerm += T(adjointRAS->FISensitivityTerm());

    // Term 7, term from objective functions
    PtrList<objective>& functions = objectiveManager_.getObjectiveFunctions();

    for (objective& objI : functions)
    {
        if (objI.hasGradDxDbMult())
        {
            flowTerm += objI.weight()*objI.gradDxDbMultiplier();
        }
    }

    flowTerm.correctBoundaryConditions();

    gradDxDbMult += flowTerm.T()*dt;
  //profiling::writeNow();
}


void Foam::incompressibleAdjointSolver::accumulateDivDxDbMultiplier
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


void Foam::incompressibleAdjointSolver::accumulateGeometryVariationsMultipliers
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
            const scalar wei = func.weight();
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


void Foam::incompressibleAdjointSolver::accumulateBCSensitivityIntegrand
(
    autoPtr<boundaryVectorField>& bcDxDbMult,
    const labelHashSet& sensitivityPatchIDs,
    const scalar dt
)
{
    if (!hasBCdxdbMult(sensitivityPatchIDs))
    {
        return;
    }

    // Grab references
    const volVectorField& Ua = getAdjointVars().Ua();
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        getAdjointVars().adjointTurbulence();

    // Fields needed to calculate adjoint sensitivities
    const autoPtr<incompressible::RASModelVariables>&
       turbVars = primalVars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = primalVars_.laminarTransport();
    volScalarField nuEff(lamTransp.nu() + turbVars->nut());
    tmp<volTensorField> tgradUa = fvc::grad(Ua);
    const volTensorField::Boundary& gradUabf = tgradUa.cref().boundaryField();

    // Avoid updating the event number to keep consistency with cases caching
    // gradUa
    auto& UaBoundary = getAdjointVars().Ua().boundaryFieldRef(false);
    auto& nuEffBoundary = nuEff.boundaryField();

    for (const label patchI : sensitivityPatchIDs)
    {
        fvPatchVectorField& Uab = UaBoundary[patchI];
        if (isA<adjointVectorBoundaryCondition>(Uab))
        {
            tmp<tensorField> dxdbMult =
                refCast<adjointVectorBoundaryCondition>(Uab).dxdbMult();
            if (dxdbMult)
            {
                const fvPatch& patch = mesh_.boundary()[patchI];
                tmp<vectorField> tnf = patch.nf();
                const scalarField& magSf = patch.magSf();

                tmp<vectorField> DvDbMult =
                    nuEffBoundary[patchI]
                   *(Uab.snGrad() + (gradUabf[patchI] & tnf))
        //        - (nf*pa.boundaryField()[patchI])
                  + adjointTurbulence().adjointMomentumBCSource()[patchI];
                bcDxDbMult()[patchI] += (DvDbMult & dxdbMult())*magSf*dt;
            }
        }
    }
}


void Foam::incompressibleAdjointSolver::accumulateOptionsDxDbMultiplier
(
    vectorField& optionsDxDbMult,
    const scalar dt
)
{
    // Terms from fvOptions - missing contributions from turbulence models
    const incompressibleAdjointVars& av = getAdjointVars();
    vectorField temp(mesh_.nCells(), Zero);
    fv::options::New(this->mesh_).postProcessSens
    (
        temp, av.UaInst().name(), av.solverName()
    );
    optionsDxDbMult += temp*dt;
    temp = Zero;
    fv::options::New(this->mesh_).postProcessSens
    (
        temp, av.paInst().name(), av.solverName()
    );
    optionsDxDbMult += temp*dt;
}


void Foam::incompressibleAdjointSolver::topOSensMultiplier
(
    scalarField& betaMult,
    const word& designVariablesName,
    const scalar dt
)
{
    const incompressibleAdjointVars& adjointVars = getAdjointVars();
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars.Ua();
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars.adjointTurbulence();
    fv::options& fvOptions(fv::options::New(mesh_));

    // Term from the momentum equations
    scalarField momSens((U.primitiveField() & Ua.primitiveField())*dt);
    Foam::sensitivityTopO::postProcessSens
        (betaMult, momSens, fvOptions, U.name(), designVariablesName);
    if (debug > 2)
    {
        volScalarField IvSens
        (
            IOobject
            (
                "IvSens" + solverName(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        IvSens.primitiveFieldRef() = momSens;
        IvSens.write();
    }

    // Term from the turbulence model.
    // Includes already the derivative of the interpolation function
    betaMult +=
        (adjointTurbulence->topologySensitivities(designVariablesName))*dt;

    // Terms resulting directly from the objective function
    PtrList<objective>& functions = objectiveManager_.getObjectiveFunctions();
    for (objective& objI : functions)
    {
        const scalar weight(objI.weight());
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
