/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "adjointSimple.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointSimple, 0);
    addToRunTimeSelectionTable
    (
        incompressibleAdjointSolver,
        adjointSimple,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::incompressibleAdjointVars& Foam::adjointSimple::allocateVars()
{
    vars_.reset
    (
        new incompressibleAdjointVars
        (
            mesh_,
            solverControl_(),
            objectiveManager_,
            primalVars_
        )
    );
    return getAdjointVars();
}


void Foam::adjointSimple::continuityErrors()
{
    const surfaceScalarField& phia = adjointVars_.phiaInst();
    volScalarField contErr(fvc::div(phia));

    scalar sumLocalContErr = mesh_.time().deltaTValue()*
        mag(contErr)().weightedAverage(mesh_.V()).value();

    scalar globalContErr = mesh_.time().deltaTValue()*
        contErr.weightedAverage(mesh_.V()).value();
    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::adjointSimple::preCalculateSensitivities()
{
    adjointSensitivity_->accumulateIntegrand(scalar(1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointSimple::adjointSimple
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
:
    incompressibleAdjointSolver
    (
        mesh,
        managerType,
        dict,
        primalSolverName,
        solverName
    ),
    solverControl_(SIMPLEControl::New(mesh, managerType, *this)),
    adjointVars_(allocateVars()),
    cumulativeContErr_(Zero)
{
    ATCModel_.reset
    (
        ATCModel::New
        (
            mesh,
            primalVars_,
            adjointVars_,
            dict.subDict("ATCModel")
        ).ptr()
    );

    setRefCell
    (
        adjointVars_.paInst(),
        solverControl_().dict(),
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
    allocateSensitivities();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointSimple::solveIter()
{
    solverControl_().incrementIter();
    if (solverControl_().performIter())
    {
        preIter();
        mainIter();
        postIter();
    }
}


void Foam::adjointSimple::preIter()
{
    Info<< "Time = " << mesh_.time().timeName() << "\n" << endl;
}


void Foam::adjointSimple::mainIter()
{
    addProfiling(adjointSimple, "adjointSimple::mainIter");
    // Grab primal references
    const surfaceScalarField& phi = primalVars_.phi();
    // Grab adjoint references
    volScalarField& pa = adjointVars_.paInst();
    volVectorField& Ua = adjointVars_.UaInst();
    surfaceScalarField& phia = adjointVars_.phiaInst();
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();
    const label&  paRefCell  = solverControl_().pRefCell();
    const scalar& paRefValue = solverControl_().pRefValue();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Momentum predictor
    //~~~~~~~~~~~~~~~~~~~

    tmp<fvVectorMatrix> tUaEqn
    (
        fvm::div(-phi, Ua)
      + adjointTurbulence->divDevReff(Ua)
      + adjointTurbulence->adjointMeanFlowSource()
      ==
        fvOptions(Ua)
    );
    fvVectorMatrix& UaEqn = tUaEqn.ref();

    // Add sources from boundary conditions
    UaEqn.boundaryManipulate(Ua.boundaryFieldRef());

    // Add sources from volume-based objectives
    objectiveManager_.addSource(UaEqn);

    // Add ATC term
    ATCModel_->addATC(UaEqn);

    // Additional source terms (e.g. energy equation)
    addMomentumSource(UaEqn);

    UaEqn.relax();

    fvOptions.constrain(UaEqn);

    if (solverControl_().momentumPredictor())
    {
        Foam::solve(UaEqn == -fvc::grad(pa));

        fvOptions.correct(Ua);
    }

    // Pressure Eq
    //~~~~~~~~~~~~
    {
        volScalarField rAUa(1.0/UaEqn.A());
        // 190402: Vag: to be updated.
        // Probably a constrainHabyA by class is needed?
        volVectorField HabyA(constrainHbyA(rAUa*UaEqn.H(), Ua, pa));
        surfaceScalarField phiaHbyA("phiaHbyA", fvc::flux(HabyA));
        adjustPhi(phiaHbyA, Ua, pa);

        tmp<volScalarField> rAtUa(rAUa);

        if (solverControl_().consistent())
        {
            rAtUa = 1.0/(1.0/rAUa - UaEqn.H1());
            phiaHbyA +=
                fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh_.magSf();
            HabyA -= (rAUa - rAtUa())*fvc::grad(pa);
        }

        tUaEqn.clear();

        // Update the pressure BCs to ensure flux consistency
        // constrainPressure(p, U, phiHbyA, rAtU(), MRF_);

        // Non-orthogonal pressure corrector loop
        while (solverControl_().correctNonOrthogonal())
        {
            fvScalarMatrix paEqn
            (
                fvm::laplacian(rAtUa(), pa) == fvc::div(phiaHbyA)
            );

            paEqn.boundaryManipulate(pa.boundaryFieldRef());

            addPressureSource(paEqn);

            fvOptions.constrain(paEqn);
            paEqn.setReference(paRefCell, paRefValue);

            paEqn.solve();

            if (solverControl_().finalNonOrthogonalIter())
            {
                phia = phiaHbyA - paEqn.flux();
            }
        }

        continuityErrors();

        // Explicitly relax pressure for adjoint momentum corrector
        pa.relax();

        // Momentum corrector
        Ua = HabyA - rAtUa()*fvc::grad(pa);
        Ua.correctBoundaryConditions();
        fvOptions.correct(Ua);
        pa.correctBoundaryConditions();
    }

    adjointTurbulence->correct();

    if (solverControl_().printMaxMags())
    {
        dimensionedScalar maxUa = gMax(mag(Ua)());
        dimensionedScalar maxpa = gMax(mag(pa)());
        Info<< "Max mag (" << Ua.name() << ") = " << maxUa.value() << endl;
        Info<< "Max mag (" << pa.name() << ") = " << maxpa.value() << endl;
    }
}


void Foam::adjointSimple::postIter()
{
    solverControl_().write();

    // Average fields if necessary
    adjointVars_.computeMeanFields();

    // Print execution time
    mesh_.time().printExecutionTime(Info);
}


void Foam::adjointSimple::solve()
{
    addProfiling(adjointSimple, "adjointSimple::solve");
    if (active_)
    {
        preLoop();
        while (solverControl_().loop())
        {
            solveIter();
        }
        postLoop();
    }
}


bool Foam::adjointSimple::loop()
{
    return solverControl_().loop();
}


void Foam::adjointSimple::preLoop()
{
    // Reset initial and mean fields before solving
    adjointVars_.restoreInitValues();
    adjointVars_.resetMeanFields();
}


void Foam::adjointSimple::addMomentumSource(fvVectorMatrix& matrix)
{
    // Does nothing
}


void Foam::adjointSimple::addPressureSource(fvScalarMatrix& matrix)
{
    // Does nothing
}


void Foam::adjointSimple::updatePrimalBasedQuantities()
{
    incompressibleAdjointSolver::updatePrimalBasedQuantities();

    // Update objective function related quantities
    objectiveManager_.updateAndWrite();
}


void Foam::adjointSimple::addTopOFvOptions() const
{
    // Determine number of variables related to the adjoint turbulence model
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();
    const wordList& turbVarNames =
        adjointTurbulence().getAdjointTMVariablesBaseNames();
    label nTurbVars(turbVarNames.size());
    if (adjointTurbulence().includeDistance())
    {
        nTurbVars++;
    }

    // Determine names of fields to be added to the dictionary
    wordList names(1 + nTurbVars);
    label varID(0);
    names[varID++] = adjointVars_.UaInst().name();
    for (const word& turbName : turbVarNames)
    {
        names[varID++] = turbName;
    }
    if (adjointTurbulence().includeDistance())
    {
        names[varID++] =
            word(useSolverNameForFields() ? "da" + solverName_ : "da");
    }

    // Add entries to dictionary
    const word dictName("topOSource" + solverName_);
    dictionary optionDict(dictName);
    optionDict.add<word>("type", "topOSource");
    optionDict.add<wordList>("names", names);
    optionDict.add<word>("function", "linear");
    optionDict.add<word>("interpolationField", "beta");

    // Construct and append fvOption
    fv::optionList& fvOptions(fv::options::New(this->mesh_));
    fvOptions.push_back(fv::option::New(dictName, optionDict, mesh_));
}


bool Foam::adjointSimple::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());
    return adjointSolver::writeData(os);
}


// ************************************************************************* //
