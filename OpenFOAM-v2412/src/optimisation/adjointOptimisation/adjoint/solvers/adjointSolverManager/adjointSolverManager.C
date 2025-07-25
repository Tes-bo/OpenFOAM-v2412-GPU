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

#include "adjointSolverManager.H"
#include "primalSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointSolverManager, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointSolverManager::adjointSolverManager
(
    fvMesh& mesh,
    autoPtr<designVariables>& designVars,
    const word& managerType,
    const dictionary& dict,
    bool overrideUseSolverName
)
:
    regIOobject
    (
        IOobject
        (
            "adjointSolverManager" + dict.dictName(),
            mesh.time().system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    mesh_(mesh),
    dict_(dict),
    managerName_(dict.dictName()),
    managerType_(managerType),
    primalSolverName_(dict.get<word>("primalSolver")),
    adjointSolvers_(0),
    objectiveSolverIDs_(0),
    oneSidedConstraintSolverIDs_(0),
    doubleSidedConstraintSolverIDs_(0),
    operatingPointWeight_
    (
        dict.getOrDefault<scalar>("operatingPointWeight", 1)
    ),
    nActiveAdjointSolvers_(0),
    designVars_(designVars)
{
    dictionary& adjointSolversDict =
        const_cast<dictionary&>(dict.subDict("adjointSolvers"));

    const wordList adjSolverNames = adjointSolversDict.toc();
    adjointSolvers_.setSize(adjSolverNames.size());
    objectiveSolverIDs_.setSize(adjSolverNames.size());
    oneSidedConstraintSolverIDs_.setSize(adjSolverNames.size());
    doubleSidedConstraintSolverIDs_.setSize(adjSolverNames.size());
    label nObjectives(0);
    label nOneSidedConstraints(0);
    label nDoubleSidedConstraints(0);
    forAll(adjSolverNames, namei)
    {
        dictionary& solverDict =
            adjointSolversDict.subDict(adjSolverNames[namei]);
        if (overrideUseSolverName)
        {
            solverDict.add<bool>("useSolverNameForFields", true);
        }
        adjointSolvers_.set
        (
            namei,
            adjointSolver::New
            (
                mesh_,
                managerType,
                solverDict,
                primalSolverName_,
                adjSolverNames[namei]
            )
        );
        if (adjointSolvers_[namei].active())
        {
            nActiveAdjointSolvers_++;
        }
        if (adjointSolvers_[namei].isDoubleSidedConstraint())
        {
            doubleSidedConstraintSolverIDs_[nDoubleSidedConstraints++] = namei;
        }
        else if (adjointSolvers_[namei].isConstraint())
        {
            oneSidedConstraintSolverIDs_[nOneSidedConstraints++] = namei;
        }
        else
        {
            objectiveSolverIDs_[nObjectives++] = namei;
        }
    }
    objectiveSolverIDs_.setSize(nObjectives);
    oneSidedConstraintSolverIDs_.setSize(nOneSidedConstraints);
    doubleSidedConstraintSolverIDs_.setSize(nDoubleSidedConstraints);

    Info<< "Found " << nOneSidedConstraints
        << " adjoint solvers acting as single-sided constraints" << endl;

    Info<< "Found " << nDoubleSidedConstraints
        << " adjoint solvers acting as double-sided constraints" << endl;

    Info<< "Found " << nActiveAdjointSolvers_
        << " active adjoint solvers" << endl;

    // Having more than one non-aggregated objectives per operating point
    // is needlessly expensive. Issue a warning
    if (objectiveSolverIDs_.size() > 1)
    {
        WarningInFunction
            << "Number of adjoint solvers corresponding to objectives "
            << "is greater than 1 (" << objectiveSolverIDs_.size() << ")" << nl
            << "Consider aggregating your objectives to one\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adjointSolverManager::readDict(const dictionary& dict)
{
    dict_ = dict;

    const dictionary& adjointSolversDict = dict.subDict("adjointSolvers");

    // Note: only updating existing solvers
    for (adjointSolver& solver : adjointSolvers_)
    {
        solver.readDict(adjointSolversDict.subDict(solver.name()));
    }

    return true;
}


const Foam::word& Foam::adjointSolverManager::managerName() const
{
    return managerName_;
}


const Foam::word& Foam::adjointSolverManager::primalSolverName() const
{
    return primalSolverName_;
}


const Foam::dictionary& Foam::adjointSolverManager::dict() const
{
    return dict_;
}


const Foam::PtrList<Foam::adjointSolver>&
Foam::adjointSolverManager::adjointSolvers() const
{
    return adjointSolvers_;
}


Foam::PtrList<Foam::adjointSolver>&
Foam::adjointSolverManager::adjointSolvers()
{
    return adjointSolvers_;
}


Foam::wordList Foam::adjointSolverManager::adjointSolversNames() const
{
    wordList names(adjointSolvers_.size());
    forAll(adjointSolvers_, sI)
    {
        names[sI]  = adjointSolvers_[sI].name();
    }
    return names;
}


Foam::scalar Foam::adjointSolverManager::operatingPointWeight() const
{
    return operatingPointWeight_;
}


Foam::label Foam::adjointSolverManager::nActiveAdjointSolvers() const
{
    return nActiveAdjointSolvers_;
}


Foam::label Foam::adjointSolverManager::nActiveAdjointSolvers
(
    const dictionary& dict
)
{
    const dictionary& adjointSolversDict = dict.subDict("adjointSolvers");
    const wordList adjSolverNames = adjointSolversDict.toc();
    label n(0);
    Switch active(true);
    forAll(adjSolverNames, namei)
    {
        active = adjointSolversDict.subDict(adjSolverNames[namei]).
            getOrDefault<bool>("active", true);
        if (active)
        {
            n++;
        }
    }
    return n;
}


Foam::label Foam::adjointSolverManager::nConstraints() const
{
    return nOneSidedConstraints() + 2*nDoubleSidedConstraints();
}


Foam::label Foam::adjointSolverManager::nOneSidedConstraints() const
{
    return oneSidedConstraintSolverIDs_.size();
}


Foam::label Foam::adjointSolverManager::nDoubleSidedConstraints() const
{
    return doubleSidedConstraintSolverIDs_.size();
}


Foam::label Foam::adjointSolverManager::nObjectives() const
{
    return objectiveSolverIDs_.size();
}


Foam::label Foam::adjointSolverManager::nAdjointSolvers() const
{
    return nOneSidedConstraints() + nDoubleSidedConstraints() + nObjectives();
}


void Foam::adjointSolverManager::solveAdjointEquations()
{
    //  Solve all adjoint equations of this adjointSolverManager
    for (adjointSolver& solver : adjointSolvers_)
    {
        // Update all primal-based quantities needed by the adjoint equations
        solver.updatePrimalBasedQuantities();

        // Solve the adjoint equations taking into consideration the weighted
        // contribution of possibly multiple objectives
        solver.solve();

        // Compute sensitivities and force writing to the adjoint dictionary
        // if this an output time
        solver.computeObjectiveSensitivities(designVars_);
        if (mesh_.time().writeTime())
        {
            solver.regIOobject::write(true);
        }
    }
}


Foam::tmp<Foam::scalarField>
Foam::adjointSolverManager::aggregateSensitivities()
{
    auto tsens = tmp<scalarField>::New();
    auto& sens = tsens.ref();

    // Sum sensitivities from all objectives expect the constraints
    for (const label solveri : objectiveSolverIDs_)
    {
        // Sum contributions
        const scalarField& solverSens =
            adjointSolvers_[solveri].getObjectiveSensitivities(designVars_);

        if (sens.empty())
        {
            sens = scalarField(solverSens.size(), Zero);
        }
        sens += solverSens;
    }

    return tsens;
}


Foam::PtrList<Foam::scalarField>
Foam::adjointSolverManager::constraintSensitivities()
{
    PtrList<scalarField> constraintSens(nConstraints());
    // Only one-sided constraints
    label cI(0);
    for (const label consI : oneSidedConstraintSolverIDs_)
    {
        constraintSens.set
        (
            cI++,
            new scalarField
                (adjointSolvers_[consI].getObjectiveSensitivities(designVars_))
        );
    }

    // Two-sided constraints. Negated left-most side sensitivities
    for (const label consI : doubleSidedConstraintSolverIDs_)
    {
        scalarField sens
            (adjointSolvers_[consI].getObjectiveSensitivities(designVars_));
        constraintSens.set(cI++, new scalarField(  sens));
        constraintSens.set(cI++, new scalarField(- sens));
    }

    return constraintSens;
}


void Foam::adjointSolverManager::computeAllSensitivities()
{
    for (adjointSolver& adjSolver : adjointSolvers_)
    {
        adjSolver.computeObjectiveSensitivities(designVars_);
    }
}


void Foam::adjointSolverManager::clearSensitivities()
{
    for (adjointSolver& adjSolver : adjointSolvers_)
    {
        adjSolver.clearSensitivities();
    }
}


Foam::scalar Foam::adjointSolverManager::objectiveValue()
{
    scalar objValue(Zero);
    for (const label solveri : objectiveSolverIDs_)
    {
        objectiveManager& objManager =
            adjointSolvers_[solveri].getObjectiveManager();
        objValue += objManager.print();
    }

    return objValue;
}


Foam::tmp<Foam::scalarField> Foam::adjointSolverManager::constraintValues()
{
    auto tconstraintValues(tmp<scalarField>::New(nConstraints(), Zero));
    scalarField& constraintValues = tconstraintValues.ref();
    label cI(0);
    // One-sided constraints only
    for (const label consI : oneSidedConstraintSolverIDs_)
    {
        objectiveManager& objManager =
            adjointSolvers_[consI].getObjectiveManager();
        constraintValues[cI++] = objManager.print();
    }
    // Double-sided constraints
    // Objective value of the left-most side is negated
    for (const label consI : doubleSidedConstraintSolverIDs_)
    {
        objectiveManager& objManager =
            adjointSolvers_[consI].getObjectiveManager();
        constraintValues[cI++] = objManager.print(false);
        constraintValues[cI++] = objManager.print(true);
    }

    return tconstraintValues;
}


void Foam::adjointSolverManager::updatePrimalBasedQuantities(const word& name)
{
    if (primalSolverName_ == name)
    {
        for (adjointSolver& solver : adjointSolvers_)
        {
            solver.updatePrimalBasedQuantities();
        }
    }
}


bool Foam::adjointSolverManager::isMaster() const
{
    return mesh_.lookupObject<primalSolver>(primalSolverName_).isMaster();
}


// ************************************************************************* //
