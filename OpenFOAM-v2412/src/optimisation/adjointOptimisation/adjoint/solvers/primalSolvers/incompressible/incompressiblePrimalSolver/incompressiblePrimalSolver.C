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

#include "incompressiblePrimalSolver.H"
#include "adjustPhi.H"
#include "adjointSolver.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressiblePrimalSolver, 0);
    defineRunTimeSelectionTable(incompressiblePrimalSolver, dictionary);
    addToRunTimeSelectionTable
    (
        primalSolver,
        incompressiblePrimalSolver,
        primalSolver
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressiblePrimalSolver::incompressiblePrimalSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& solverName
)
:
    primalSolver(mesh, managerType, dict, solverName),
    phiReconstructionTol_
    (
        dict.subOrEmptyDict("fieldReconstruction").
            getOrDefault<scalar>("tolerance", 5.e-5)
    ),
    phiReconstructionIters_
    (
        dict.subOrEmptyDict("fieldReconstruction").
            getOrDefault<label>("iters", 10)
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::incompressiblePrimalSolver>
Foam::incompressiblePrimalSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
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
            "incompressiblePrimalSolver",
            solverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return
        autoPtr<incompressiblePrimalSolver>
        (
            ctorPtr(mesh, managerType, dict, solverName)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::incompressiblePrimalSolver::readDict(const dictionary& dict)
{
    if (primalSolver::readDict(dict))
    {
        return true;
    }

    return false;
}


Foam::UPtrList<Foam::objective>
Foam::incompressiblePrimalSolver::getObjectiveFunctions() const
{
    DynamicList<objective*> objectives;

    for (auto& adjoint : mesh_.sorted<adjointSolver>())
    {
        if (adjoint.primalSolverName() == solverName_)
        {
            PtrList<objective>& managerObjectives =
                adjoint.getObjectiveManager().getObjectiveFunctions();

            for (objective& obj : managerObjectives)
            {
                objectives.push_back(&obj);
            }
        }
    }

    return UPtrList<objective>(objectives);
}


const Foam::incompressibleVars&
Foam::incompressiblePrimalSolver::getIncoVars() const
{
    const incompressibleVars& incoVars =
        refCast<incompressibleVars>(const_cast<variablesSet&>(vars_()));
    return incoVars;
}


Foam::incompressibleVars&
Foam::incompressiblePrimalSolver::getIncoVars()
{
    incompressibleVars& incoVars =
        refCast<incompressibleVars>(const_cast<variablesSet&>(vars_()));
    return incoVars;
}


void Foam::incompressiblePrimalSolver::correctBoundaryConditions()
{
    incompressibleVars& vars = getIncoVars();
    // Update boundary conditions for all primal volFields,
    // including averaged ones, if present
    vars.correctBoundaryConditions();

    // phi cannot be updated through correctBoundaryConditions.
    // Re-compute based on the Rhie-Chow interpolation scheme.
    // This is a non-linear process
    // (phi depends on UEqn().A() which depends on phi)
    // so iterations are required

    volScalarField& p = vars.p();
    volVectorField& U = vars.U();
    surfaceScalarField& phi = vars.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence = vars.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    scalar contError(GREAT), diff(GREAT);
    for (label iter = 0; iter < phiReconstructionIters_; ++iter)
    {
        Info<< "phi correction iteration " << iter << endl;

        // Form momentum equations to grab A
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          + turbulence->divDevReff(U)
          ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();
        UEqn.relax();
        fvOptions.constrain(UEqn);

        // Pressure equation will give the Rhie-Chow correction
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA("HbyA", U);
        HbyA = rAU*UEqn.H();
        tUEqn.clear();

        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        adjustPhi(phiHbyA, U, p);

        //fvOptions.makeRelative(phiHbyA);

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );
        phi = phiHbyA - pEqn.flux();

        // Check convergence
        // ~~~~~~~~~~~~~~~~~
        scalar contErrorNew =
            mesh_.time().deltaTValue()*
            mag(fvc::div(phi)())().weightedAverage(mesh_.V()).value();

        Info<< "sum local = " << contErrorNew << endl;
        diff = mag(contErrorNew - contError)/contError;
        contError = contErrorNew;

        if (diff < phiReconstructionTol_) break;

        Info<< endl;
     }
}


// ************************************************************************* //
