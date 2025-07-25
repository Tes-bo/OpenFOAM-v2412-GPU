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

#include "adjointSolverManager.H"
#include "dictionary.H"
#include "optimisationManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(optimisationManager, 0);
    defineRunTimeSelectionTable(optimisationManager, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::optimisationManager::resetTime
(
    const dimensionedScalar startTime,
    const label startTimeIndex,
    const scalar endTime
)
{
    // Does nothing in base
}


void Foam::optimisationManager::lineSearchUpdate()
{
    // Compute direction of update
    tmp<scalarField> tdirection = dvUpdate_->computeDirection();
    const scalarField& direction = tdirection.ref();

    // Grab reference to line search
    autoPtr<lineSearch>& lineSrch = dvUpdate_->getLineSearch();

    // Store starting point
    dvUpdate_->getDesignVariables()->storeDesignVariables();

    // Compute merit function before update
    scalar meritFunction = dvUpdate_->computeMeritFunction();
    lineSrch->setOldMeritValue(meritFunction);
    dvUpdate_->setOldObjectiveValue();

    // Get merit function derivative
    scalar dirDerivative =
        dvUpdate_->meritFunctionDirectionalDerivative();
    lineSrch->setDeriv(dirDerivative);
    lineSrch->setDirection(direction);

    // Reset initial step.
    // Might be interpolated from previous optimisation cycles
    lineSrch->reset();

    // Perform line search
    while (lineSrch->loop())
    {
        Info<< "\n- - - - - - - - - - - - - - -"  << endl;
        Info<< "Line search iteration "   << lineSrch->innerIter() << endl;
        Info<< "- - - - - - - - - - - - - - -\n"  << endl;

        // Update design variables. Multiplication with line search step
        // happens inside the update(direction) function
        moveDesignVariables(direction);

        const dimensionedScalar startTime = mesh_.time();
        const label startTimeIndex = mesh_.time().timeIndex();
        const scalar primalEndTime = mesh_.time().endTime().value();
        // Solve all primal equations
        solvePrimalEquations();

        // Compute and set new merit function
        meritFunction = dvUpdate_->computeMeritFunction();
        lineSrch->setNewMeritValue(meritFunction);

        if (lineSrch->computeGradient())
        {
            // Reset adjoint sensitivities in all adjoint solver managers
            clearSensitivities();

            // Solve all adjoint equations
            solveAdjointEquations();

            // Update objective and gradient information known by updateMethod
            dvUpdate_->updateGradientsAndValues();

            // Update the directional derivative
            dirDerivative =
                dvUpdate_->meritFunctionDirectionalDerivative();
            lineSrch->setNewDeriv(dirDerivative);
        }

        if (lineSrch->converged())
        {
            // If line search criteria have been met, proceed
            Info<< "Line search converged in " << lineSrch->innerIter()
                << " iterations." << endl;
            scalarField scaledCorrection(lineSrch->step()*direction);
            dvUpdate_->postUpdate(scaledCorrection);
            break;
        }
        else
        {
            // If maximum number of iteration has been reached, continue
            if (lineSrch->innerIter() == lineSrch->maxIters())
            {
                Info<< "Line search reached max. number of iterations.\n"
                    << "Proceeding to the next optimisation cycle" << endl;
                scalarField scaledCorrection(lineSrch->step()*direction);
                dvUpdate_->postUpdate(scaledCorrection);
            }
            // Reset to initial design variables and update step
            else
            {
                //- Reset time if necessary
                this->resetTime(startTime, startTimeIndex, primalEndTime);
                dvUpdate_->getDesignVariables()->resetDesignVariables();
                lineSrch->updateStep();
            }
        }
    }

    // If line search did not need to get the new gradient, do it now in order
    // to have it for the next optimisation cycle
    if (!lineSrch->computeGradient())
    {
        // Reset adjoint sensitivities in all adjoint solver managers
        clearSensitivities();

        // Solve all adjoint equations
        solveAdjointEquations();
    }
}


void Foam::optimisationManager::fixedStepUpdate()
{
    // Update design variables
    if (shouldUpdateDesignVariables_)
    {
        moveDesignVariables();
    }

    // Solve primal equations
    solvePrimalEquations();

    // Check the convergence criteria of the optimisation loop
    dvUpdate_->checkConvergence();

    // Reset adjoint sensitivities in all adjoint solver managers
    clearSensitivities();

    // Solve all adjoint equations
    solveAdjointEquations();
}


void Foam::optimisationManager::moveDesignVariables()
{
    // Update design variables
    dvUpdate_->update();
}


void Foam::optimisationManager::moveDesignVariables
(
    const scalarField& direction
)
{
    // Update design variables
    dvUpdate_->update(direction);
}


void Foam::optimisationManager::initialize()
{
    dictionary& primalSolversDict = subDict("primalSolvers");
    const wordList& primalSolverNames = primalSolversDict.toc();

    // Construct primal solvers
    primalSolvers_.setSize(primalSolverNames.size());
    forAll(primalSolvers_, solveri)
    {
        dictionary& solverDict =
            primalSolversDict.subDict(primalSolverNames[solveri]);
        if (primalSolvers_.size() > 1)
        {
            solverDict.add<bool>("useSolverNameForFields", true);
        }
        primalSolvers_.set
        (
            solveri,
            primalSolver::New
            (
                mesh_,
                managerType_,
                solverDict,
                primalSolverNames[solveri]
            )
        );
    }

    // Construct adjointSolverManagers
    const dictionary& adjointManagersDict = subDict("adjointManagers");
    const wordList adjointManagerNames = adjointManagersDict.toc();
    adjointSolverManagers_.setSize(adjointManagerNames.size());

    // Determine the number of adjoint solvers which are not null (i.e. need to
    // allocate adjoint fields). Used to determine whether the adjoint field
    // names should be appended by the solver name
    label nNotNullAdjointSolvers(0);
    for (const word& adjManager : adjointManagerNames)
    {
        const dictionary& adjSolversDict =
            adjointManagersDict.subDict(adjManager).subDict("adjointSolvers");
        const wordList adjointSolverNames = adjSolversDict.toc();
        for (const word& adjSolver : adjointSolverNames)
        {
            if (adjSolversDict.subDict(adjSolver).get<word>("type") != "null")
            {
                ++nNotNullAdjointSolvers;
            }
        }
    }
    Info<< "Found "
        << nNotNullAdjointSolvers
        << " adjoint solvers that allocate fields"
        << endl;
    bool overrideUseSolverName(nNotNullAdjointSolvers > 1);

    forAll(adjointSolverManagers_, manageri)
    {
        adjointSolverManagers_.set
        (
            manageri,
            new adjointSolverManager
            (
                mesh_,
                designVars_,
                managerType_,
                adjointManagersDict.subDict(adjointManagerNames[manageri]),
                overrideUseSolverName
            )
        );
    }

    // Sanity checks on the naming convention
    if (primalSolvers_.size() > 1)
    {
        for (const primalSolver& solveri : primalSolvers_)
        {
            if (!solveri.useSolverNameForFields())
            {
                FatalErrorInFunction
                    << "Multiple primal solvers are present but "
                    << "useSolverNameForFields is set to false in "
                    << "primal solver " << solveri.solverName() << nl
                    << "This is considered fatal."
                    << exit(FatalError);
            }
        }
    }

    if (nNotNullAdjointSolvers > 1)
    {
        for (const adjointSolverManager& amI : adjointSolverManagers_)
        {
            const PtrList<adjointSolver>& adjointSolvers = amI.adjointSolvers();
            for (const adjointSolver& asI : adjointSolvers)
            {
                if (!asI.useSolverNameForFields())
                {
                    FatalErrorInFunction
                        << "Multiple adjoint solvers are present but "
                        << "useSolverNameForFields is set to false in "
                        << "adjoint solver " << asI.solverName() << nl
                        << "This is considered fatal."
                        << exit(FatalError);
                }
            }
        }
    }
    if (designVars_)
    {
        designVars_().addFvOptions(primalSolvers_, adjointSolverManagers_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optimisationManager::optimisationManager(fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "optimisationDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    mesh_(mesh),
    time_(const_cast<Time&>(mesh.time())),
    designVars_(nullptr),
    primalSolvers_(),
    adjointSolverManagers_(),
    managerType_(get<word>("optimisationManager")),
    dvUpdate_(nullptr),
    shouldUpdateDesignVariables_(true)
{
    // The "designVariables" sub-dictionary is optional
    const dictionary odict(this->subOrEmptyDict("optimisation"));
    const dictionary* designVarsDictPtr = odict.findDict("designVariables");

    if (designVarsDictPtr)
    {
        designVars_ = designVariables::New(mesh_, *designVarsDictPtr);
    }
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::optimisationManager> Foam::optimisationManager::New
(
    fvMesh& mesh
)
{
    const IOdictionary dict
    (
        IOobject
        (
            "optimisationDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const word modelType(dict.get<word>("optimisationManager"));

    Info<< "optimisationManager type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "optimisationManager",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<optimisationManager>(ctorPtr(mesh));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool Foam::optimisationManager::read()
{
    if (regIOobject::read())
    {
        // Note: Only changing existing solvers - not adding any new
        const dictionary& primalSolversDict = subDict("primalSolvers");
        for (primalSolver& sol : primalSolvers_)
        {
            sol.readDict(primalSolversDict.subDict(sol.solverName()));
        }

        const dictionary& adjointManagersDict = subDict("adjointManagers");
        for (adjointSolverManager& man : adjointSolverManagers_)
        {
            man.readDict(adjointManagersDict.subDict(man.managerName()));
        }

        if (designVars_)
        {
            designVars_->readDict
                (subDict("optimisation").subDict("designVariables"));
        }

        return true;
    }

    return false;
}


void Foam::optimisationManager::updateDesignVariables()
{
    // Update design variables using either a line-search scheme or
    // a fixed-step update
    if (dvUpdate_->getLineSearch())
    {
        lineSearchUpdate();
    }
    else
    {
        fixedStepUpdate();
    }
}


void Foam::optimisationManager::solvePrimalEquations()
{
    // Solve all primal equations
    forAll(primalSolvers_, psI)
    {
        primalSolvers_[psI].solve();
    }
}


void Foam::optimisationManager::solveAdjointEquations()
{
    // Solve all adjoint solver equations
    forAll(adjointSolverManagers_, amI)
    {
        adjointSolverManagers_[amI].solveAdjointEquations();
    }
}


void Foam::optimisationManager::computeSensitivities()
{
    // Compute senstivities from all adjoint solvers
    forAll(adjointSolverManagers_, amI)
    {
        adjointSolverManagers_[amI].computeAllSensitivities();
    }
}


void Foam::optimisationManager::clearSensitivities()
{
    for (adjointSolverManager& adjSolvManager : adjointSolverManagers_)
    {
        adjSolvManager.clearSensitivities();
    }
}


void Foam::optimisationManager::updatePrimalBasedQuantities()
{
    forAll(adjointSolverManagers_, amI)
    {
        PtrList<adjointSolver>& adjointSolvers =
            adjointSolverManagers_[amI].adjointSolvers();

        forAll(adjointSolvers, asI)
        {
            adjointSolvers[asI].updatePrimalBasedQuantities();
        }
    }
}


// ************************************************************************* //
