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
#include "designVariablesUpdate.H"
#include "constrainedOptimisationMethod.H"
#include "adjointNull.H"
#include "IOmanip.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(designVariablesUpdate, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::designVariablesUpdate::nConstraints
(
    PtrList<adjointSolverManager>& adjointSolverManagers
) const
{
    // Figure out number of adjoint solvers corresponding to constraints.
    // Looks in all operating points
    label nConstraints(0);
    for (const adjointSolverManager& adjManagerI : adjointSolvManagers_)
    {
        nConstraints += adjManagerI.nConstraints();
    }
    // Add constraints that might emerge from the design variables
    tmp<scalarField> designVarsConstraints = designVars_().constraintValues();
    if (designVarsConstraints)
    {
        nConstraints += designVarsConstraints().size();
    }
    return nConstraints;
}


Foam::label Foam::designVariablesUpdate::nAdjointSolvers() const
{
    label n(0);
    for (adjointSolverManager& adjSolvManager : adjointSolvManagers_)
    {
        for (adjointSolver& solver : adjSolvManager.adjointSolvers())
        {
            if (!isA<adjointNull>(solver))
            {
                ++n;
            }
        }
    }
    return n;
}


void Foam::designVariablesUpdate::writeCPUcostHeader()
{
    unsigned int width(IOstream::defaultPrecision() + 5);
    CPUcostFile_
        << setw(width) << "#Cycle" << " "
        << setw(width) << "LineSearchIters" << " "
        << setw(width) << "CycleCPUcost" << " "
        << setw(width) << "CyclePrimalSolutions" << " "
        << setw(width) << "CycleAdjointSolutions" << " "
        << setw(width) << "TotalCPUcost" << " "
        << setw(width) << "TotalPrimalSolutions" << " "
        << setw(width) << "TotalAdjointSolutions" << endl;
}


void Foam::designVariablesUpdate::writeToCostFile(bool zeroAdjointSolns)
{
    unsigned int width(IOstream::defaultPrecision() + 5);
    label cyclePrimalSolutions(nPrimalsPerCycle_);
    label cycleAdjointSolutions(nAdjointsPerCycle_);
    label lineSearchIters(1);
    if (lineSearch_)
    {
        lineSearchIters = lineSearch_().innerIter();
        cyclePrimalSolutions *= lineSearchIters;
        if (lineSearch_().computeGradient())
        {
            cycleAdjointSolutions *= lineSearchIters;
        }
    }
    if (zeroAdjointSolns)
    {
        cycleAdjointSolutions = 0;
    }
    nPrimalSolutions_ += cyclePrimalSolutions;
    nAdjointSolutions_ += cycleAdjointSolutions;
    const scalar elapsedCpuTime = mesh_.time().elapsedCpuTime();
    const scalar cycleCost = elapsedCpuTime - CPUcost_;
    CPUcost_ = elapsedCpuTime;

    CPUcostFile_
        << setw(width) << mesh_.time().timeName() << " "
        << setw(width) << lineSearchIters << " "
        << setw(width) << cycleCost  << " "
        << setw(width) << cyclePrimalSolutions << " "
        << setw(width) << cycleAdjointSolutions << " "
        << setw(width) << CPUcost_ << " "
        << setw(width) << nPrimalSolutions_ << " "
        << setw(width) << nAdjointSolutions_ << endl;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::designVariablesUpdate::designVariablesUpdate
(
    fvMesh& mesh,
    const dictionary& dict,
    PtrList<adjointSolverManager>& adjointSolverManagers,
    autoPtr<designVariables>& designVars
)
:
    mesh_(mesh),
    dict_(dict),
    adjointSolvManagers_(adjointSolverManagers),
    designVars_(designVars),
    updateMethod_
    (
        updateMethod::New
        (
            mesh_,
            dict_.subDict("updateMethod"),
            designVars_,
            nConstraints(adjointSolverManagers)
        )
    ),
    lineSearch_
    (
        lineSearch::New
        (
            dict_.subDict("updateMethod").subOrEmptyDict("lineSearch"),
            mesh.time(),
            updateMethod_.ref()
        )
    ),
    CPUcostFile_(mesh_.time().globalPath()/"optimisation"/"CPUcost"),
    nPrimalsPerCycle_(adjointSolverManagers.size()),
    nAdjointsPerCycle_(nAdjointSolvers()),
    nPrimalSolutions_(nPrimalsPerCycle_),
    nAdjointSolutions_(nAdjointsPerCycle_),
    CPUcost_(0),
    designVarsThreshold_(nullptr),
    objectiveThreshold_(nullptr),
    convergenceCriteriaDefined_(false),
    feasibilityThreshold_
    (
        dict.subOrEmptyDict("convergence").
            getOrDefault<scalar>("feasibilityThreshold", 1.e-06)
    )
{
    dictionary convergenceDict = dict.subOrEmptyDict("convergence");
    if (convergenceDict.found("designVariables"))
    {
        designVarsThreshold_.reset
            (new scalar(convergenceDict.get<scalar>("designVariables")));
    }
    if (convergenceDict.found("objective"))
    {
        objectiveThreshold_.reset
            (new scalar(convergenceDict.get<scalar>("objective")));
    }
    convergenceCriteriaDefined_ = designVarsThreshold_ || objectiveThreshold_;
    // Check whether eta of maxInitChange are set
    if (!designVars_().isMaxInitChangeSet() && !updateMethod_().initialEtaSet())
    {
        FatalErrorInFunction
            << "Neither eta (updateMethod) or maxInitChange (designVariables) "
            << "is set."
            << exit(FatalError);
    }

    label nConstr(nConstraints(adjointSolvManagers_));
    // Sanity checks for combinations of number of constraints and
    // optimisation methods
    if
    (
        nConstr
     && !isA<constrainedOptimisationMethod>(updateMethod_())
    )
    {
        const auto& cnstrTable =
            *(constrainedOptimisationMethod::dictionaryConstructorTablePtr_);

        // Has constraints but is not a constraint optimisation method
        FatalErrorInFunction
            << "Found " << nConstr << " adjoint solvers corresponding to "
            << "constraints but the optimisation method ("
            << updateMethod_().type()
            << ") is not a constrainedOptimisationMethod." << nl
            << "Available constrainedOptimisationMethods:" << nl
            << cnstrTable.sortedToc()
            << exit(FatalError);
    }
    else if
    (
        !nConstr
     && isA<constrainedOptimisationMethod>(updateMethod_())
    )
    {
        // Does not have constraints but is a constrained optimisation method
        WarningInFunction
            << "Did not find any adjoint solvers corresponding to "
            << "constraints but the optimisation method ("
            << updateMethod_().type()
            << ") is a constrainedOptimisationMethod." << nl << nl
            << "This can cause some constraintOptimisationMethods to misbehave."
            << nl << nl
            << "Either the isConstraint bool is not set in one of the adjoint "
            << "solvers or you should consider using an updateMethod "
            << "that is not a constrainedOptimisationMethod"
            << nl << endl;
    }

    if (designVarsThreshold_)
    {
        Info<< "Optimisation will run until the max. scaled correction "
            << "of the design variables is < " << designVarsThreshold_()
            << endl;
    }
    if (objectiveThreshold_)
    {
        Info<< "Optimisation will run until the relative update of the "
            << "objective function is < " << objectiveThreshold_()
            << endl;
    }
    if (!convergenceCriteriaDefined_)
    {
        Info<< "No convergence criterion defined for optimsation" << nl
            << "It will run for " << mesh_.time().endTime().value()
            << " optimisation cycles " << nl << endl;
    }
    Info<< "Feasibility threshold is " << feasibilityThreshold_ << endl;

    // Write header of the cost file
    writeCPUcostHeader();
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::designVariablesUpdate::update()
{
    // Compute update of the design variables
    tmp<scalarField> tcorrection(computeDirection());
    scalarField& correction = tcorrection.ref();

    // Set the old value of the objective function
    setOldObjectiveValue();

    // Update the design variables
    designVars_->update(correction);

    // If direction has been scaled (say by setting the initial eta), the
    // old correction has to be updated
    postUpdate(correction);
}


void Foam::designVariablesUpdate::update(const scalarField& direction)
{
    // Multiply with line search step, if necessary
    scalarField correction(direction);
    if (lineSearch_)
    {
        lineSearch_->updateCorrection(correction);
    }

    // Update the design variables
    designVars_->update(correction);
}


Foam::tmp<Foam::scalarField> Foam::designVariablesUpdate::computeDirection()
{
    updateGradientsAndValues();
    updateMethod_->computeCorrection();
    scalarField& correction = updateMethod_->returnCorrection();

    // Compute eta if needed
    if (!updateMethod_->initialEtaSet() || designVars_->resetEta())
    {
        const scalar eta(designVars_->computeEta(correction));
        updateMethod_->modifyStep(eta);
        updateMethod_->initialEtaSet() = true;
    }

    // Intentionally copies result to new field
    return tmp<scalarField>::New(correction);
}


void Foam::designVariablesUpdate::updateGradientsAndValues()
{
    scalarField objectiveSens;
    PtrList<scalarField> constraintSens;
    scalar objectiveValue(Zero);
    DynamicList<scalar> constraintValues;

    for (adjointSolverManager& adjSolvManager : adjointSolvManagers_)
    {
        const scalar opWeight = adjSolvManager.operatingPointWeight();

        // Aggregate sensitivities of solvers corresponding to objectives
        // (i.e. not constraints)
        tmp<scalarField> tadjointSolverManagerSens =
            adjSolvManager.aggregateSensitivities();

        // Aggregate objective values of solvers corresponding to objectives
        // (i.e. not constraints)
        objectiveValue += opWeight*adjSolvManager.objectiveValue();

        // Get constraint sensitivities
        PtrList<scalarField> adjointSolverManagerConstSens =
            adjSolvManager.constraintSensitivities();

        // Post process sensitivities if needed.
        // Done here since each adjointSolverManager might post-process
        // its sensitivities in a different way
        designVars_->postProcessSens
        (
            tadjointSolverManagerSens.ref(),
            adjointSolverManagerConstSens,
            adjSolvManager.adjointSolversNames(),
            adjSolvManager.isMaster()
        );

        if (objectiveSens.empty())
        {
            objectiveSens.setSize(tadjointSolverManagerSens().size(), Zero);
        }

        // Accumulate sensitivities
        objectiveSens += opWeight*tadjointSolverManagerSens();
        forAll(adjointSolverManagerConstSens, sI)
        {
            constraintSens.
                push_back(adjointSolverManagerConstSens.set(sI, nullptr));
        }
        constraintValues.push_back(adjSolvManager.constraintValues());
    }
    // Add contraint values and gradients from the design variables
    tmp<scalarField> designVarsConstValues = designVars_->constraintValues();
    PtrList<scalarField> designVarsConstDerivs =
        designVars_->constraintDerivatives();
    if (designVarsConstValues && designVarsConstDerivs.size())
    {
        if (designVarsConstValues().size() != designVarsConstDerivs.size())
        {
            FatalErrorInFunction
                << "Size of design variables constraints and derivatives differ"
                << endl
                << exit(FatalError);
        }
        constraintValues.push_back(designVarsConstValues());
        constraintSens.push_back(std::move(designVarsConstDerivs));
    }

    // Update objective/constraint values/gradients, known by the update method
    updateMethod_->setObjectiveDeriv(objectiveSens);
    updateMethod_->setConstraintDeriv(constraintSens);
    updateMethod_->setObjectiveValue(objectiveValue);
    updateMethod_->setConstraintValues
        (scalarField(std::move(constraintValues)));
}


Foam::scalar Foam::designVariablesUpdate::computeObjectiveValue()
{
    scalar objectiveValue(Zero);
    for (adjointSolverManager& adjSolvManager : adjointSolvManagers_)
    {
        const scalar opWeight = adjSolvManager.operatingPointWeight();
        objectiveValue += opWeight*adjSolvManager.objectiveValue();
    }
    return objectiveValue;
}


void Foam::designVariablesUpdate::setOldObjectiveValue()
{
    updateMethod_->setObjectiveValueOld(computeObjectiveValue());
}


Foam::scalar Foam::designVariablesUpdate::computeMeritFunction()
{
    // Compute new objective and constraint values and update the ones
    // in updateMethod
    scalar objectiveValue(Zero);
    DynamicList<scalar> constraintValues;

    for (adjointSolverManager& adjSolvManager : adjointSolvManagers_)
    {
        const scalar opWeight = adjSolvManager.operatingPointWeight();

        objectiveValue += opWeight*adjSolvManager.objectiveValue();
        constraintValues.push_back(adjSolvManager.constraintValues());
    }

    // Add constraints directly imposed to the design variables
    tmp<scalarField> designVarsConstValues = designVars_->constraintValues();
    if (designVarsConstValues)
    {
        constraintValues.push_back(designVarsConstValues());
    }
    updateMethod_->setObjectiveValue(objectiveValue);
    updateMethod_->setConstraintValues
        (scalarField(std::move(constraintValues)));

    return updateMethod_->computeMeritFunction();
}


Foam::scalar Foam::designVariablesUpdate::meritFunctionDirectionalDerivative()
{
    return updateMethod_->meritFunctionDirectionalDerivative();
}


void Foam::designVariablesUpdate::updateOldCorrection
(
    const scalarField& oldCorrection
)
{
    updateMethod_->updateOldCorrection(oldCorrection);
}


void Foam::designVariablesUpdate::write()
{
    updateMethod_->writeAuxiliaryData();
    designVars_->writeDesignVars();
    writeToCostFile();
}


void Foam::designVariablesUpdate::postUpdate(const scalarField& oldCorrection)
{
    updateOldCorrection(oldCorrection);
    write();
    designVars_->evolveNumber();
    if (lineSearch_)
    {
        lineSearch_->postUpdate();

        // Append additional empty line at the end of the instantaneous
        // objective file to indicate the end of the block corresponding to
        // this optimisation cycle
        for (adjointSolverManager& am : adjointSolvManagers_)
        {
            for (adjointSolver& as : am.adjointSolvers())
            {
                PtrList<objective>& objectives =
                    as.getObjectiveManager().getObjectiveFunctions();
                for (objective& obj : objectives)
                {
                    obj.writeInstantaneousSeparator();
                }
            }
        }
        checkConvergence();
    }
}


void Foam::designVariablesUpdate::checkConvergence()
{
    if (!convergenceCriteriaDefined_)
    {
        return;
    }

    bool converged(false);
    const scalarField& oldCorrection = updateMethod_->returnCorrection();
    // Design variables convergence check
    if (designVarsThreshold_)
    {
        const labelList& activeVarIDs =
            designVars_->activeDesignVariables();
        const scalarField correction(oldCorrection, activeVarIDs);
        const scalarField activeVars(designVars_->getVars(), activeVarIDs);
        const scalar scaledCorrection =
            gMax(mag(correction)/(mag(activeVars) + SMALL));
        DebugInfo
            << "Current correction " << correction << nl
            << "Active vars " << activeVars << endl;
        Info<< "Max. scaled correction of the design variables = "
            << scaledCorrection
            << endl;
        if (scaledCorrection < designVarsThreshold_())
        {
            Info<< tab << "Design variables have converged " << endl;
            converged = true;
        }
    }
    // Objective convergence check
    if (objectiveThreshold_ && updateMethod_->getObjectiveValueOld())
    {
        const scalar newObjective = computeObjectiveValue();
        const scalar oldObjective = updateMethod_->getObjectiveValueOld()();
        const scalar relativeUpdate =
            mag(newObjective - oldObjective)/(mag(oldObjective) + SMALL);
        Info<< "Relative change of the objective value = "
            << relativeUpdate
            << endl;
        if (relativeUpdate < objectiveThreshold_())
        {
            Info<< tab << "Objective function has converged " << endl;
            converged = true;
        }
    }
    // Feasibility check
    const scalarField& constraints = updateMethod_->getConstraintValues();
    const scalar feasibility = sum(pos(constraints)*constraints);
    Info<< "Feasibility = " << feasibility << endl;
    if (converged && feasibility < feasibilityThreshold_)
    {
        Info<< "Stopping criteria met and all constraints satisfied." << nl
            << "Optimisation has converged, stopping ..." << nl << nl
            << "End" << nl
            << endl;
        // Force writing of all objective and constraint functions, to get
        // the converged results to files
        for (adjointSolverManager& am : adjointSolvManagers_)
        {
            for (adjointSolver& as : am.adjointSolvers())
            {
                // Use dummy weighted objective
                as.getObjectiveManager().writeObjectives();
            }
        }
        writeToCostFile(true);
        if (UPstream::parRun())
        {
            UPstream::exit(0);
        }
        else
        {
            std::exit(0);
        }
    }
}


// ************************************************************************* //
