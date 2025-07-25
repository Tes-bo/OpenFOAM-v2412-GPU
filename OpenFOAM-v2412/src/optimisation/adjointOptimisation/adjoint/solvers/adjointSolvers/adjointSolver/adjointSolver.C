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

#include "adjointSolver.H"
#include "adjointSensitivity.H"
#include "designVariables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointSolver, 0);
    defineRunTimeSelectionTable(adjointSolver, adjointSolver);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::adjointSolver::allocateSensitivities()
{
    if (computeSensitivities_)
    {
        adjointSensitivity_.reset
        (
            adjointSensitivity::New(mesh_, designVarsDict(), *this).ptr()
        );
    }
}


Foam::dictionary Foam::adjointSolver::designVarsDict() const
{
    // Re-read optimisationDict here to cover multi-region cases
    return
        IOdictionary
        (
            IOobject
            (
                "optimisationDict",
                mesh_.time().globalPath()/"system",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        ).subDict("optimisation").subDict("designVariables");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointSolver::adjointSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
:
    solver(mesh, managerType, dict, solverName),
    primalSolverName_(primalSolverName),
    objectiveManager_
    (
        mesh,
        dict.subDict("objectives"),
        solverName_,
        primalSolverName
    ),
    sensitivities_(nullptr),
    computeSensitivities_
    (
        dict.getOrDefault<bool>("computeSensitivities", true)
    ),
    isConstraint_(dict.getOrDefault<bool>("isConstraint", false)),
    isDoubleSidedConstraint_
        (dict.getOrDefault<bool>("isDoubleSidedConstraint", false)),
    adjointSensitivity_(nullptr)
{
    // Force solver to not be a (single-sided) contraint if flagged as
    // double-sided
    if (isDoubleSidedConstraint_)
    {
        isConstraint_ = false;
    }
    // Update objective-related quantities to get correct derivatives
    // in case of continuation
    objectiveManager_.update();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::adjointSolver> Foam::adjointSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
{
    const word solverType(dict.get<word>("type"));

    auto* ctorPtr = adjointSolverConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "adjointSolver",
            solverType,
            *adjointSolverConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<adjointSolver>
    (
        ctorPtr(mesh, managerType, dict, primalSolverName, solverName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adjointSolver::readDict(const dictionary& dict)
{
    if (solver::readDict(dict))
    {
        computeSensitivities_ =
            dict.getOrDefault<bool>("computeSensitivities", true);

        objectiveManager_.readDict(dict.subDict("objectives"));

        if (adjointSensitivity_)
        {
            adjointSensitivity_().readDict(designVarsDict());
        }

        return true;
    }

    return false;
}


bool Foam::adjointSolver::includeDistance() const
{
    return false;
}


Foam::dimensionSet Foam::adjointSolver::daDimensions() const
{
    NotImplemented;
    return dimless;
}


Foam::dimensionSet Foam::adjointSolver::maDimensions() const
{
    NotImplemented;
    return dimless;
}


Foam::tmp<Foam::volScalarField> Foam::adjointSolver::adjointEikonalSource()
{
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::adjointSolver::yWall() const
{
    return nullptr;
}


void Foam::adjointSolver::computeObjectiveSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    if (computeSensitivities_)
    {
        preCalculateSensitivities();
        const scalarField& sens =
            adjointSensitivity_->calculateSensitivities(designVars);
        if (!sensitivities_)
        {
            sensitivities_.reset(new scalarField(sens.size(), Zero));
        }
        sensitivities_.ref() = sens;
    }
    else
    {
        sensitivities_.reset(new scalarField());
    }
}


const Foam::scalarField& Foam::adjointSolver::getObjectiveSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    if (!sensitivities_)
    {
        // Read sensitivities from file in case of continuation
        // Done here and not in allocateSensitivities since the size of the
        // design variables and, hence, the sensitivities is not known there
        if (dictionary::found("sensitivities"))
        {
            sensitivities_ =
                tmp<scalarField>::New
                    ("sensitivities", *this, designVars().size());
        }
        else
        {
            computeObjectiveSensitivities(designVars);
        }
    }

    return sensitivities_();
}


void Foam::adjointSolver::clearSensitivities()
{
    if (computeSensitivities_)
    {
        adjointSensitivity_->clearSensitivities();
        sensitivities_.clear();
    }
}


void Foam::adjointSolver::updatePrimalBasedQuantities()
{
    // Does nothing in base
}


bool Foam::adjointSolver::writeData(Ostream& os) const
{
    if (sensitivities_.good())
    {
        sensitivities_().writeEntry("sensitivities", os);
    }
    return true;
}


// ************************************************************************* //
