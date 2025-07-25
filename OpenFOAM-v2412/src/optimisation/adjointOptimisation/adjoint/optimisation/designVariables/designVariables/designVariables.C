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

#include "designVariables.H"
#include "adjointSensitivity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(designVariables, 0);
    defineRunTimeSelectionTable(designVariables, designVariables);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::designVariables::readBounds
(
    autoPtr<scalar> lowerBoundPtr,
    autoPtr<scalar> upperBoundPtr
)
{
    // Read lower bounds for the design variables, if present
    if (dict_.found("lowerBounds"))
    {
        scalarField lowerBounds(dict_.get<scalarField>("lowerBounds"));
        if (lowerBounds.size() != getVars().size())
        {
            FatalErrorInFunction
                << "Inconsistent dimensions for lowerBounds ("
                << lowerBounds.size()
                << ") and design variables ("
                << getVars().size() << ")"
                << exit(FatalError);
        }
        lowerBounds_.reset(new scalarField(lowerBounds));
    }
    else if (dict_.found("lowerBound"))
    {
        scalar lowerBound(dict_.get<scalar>("lowerBound"));
        lowerBounds_.reset(new scalarField(getVars().size(), lowerBound));
    }
    else if (lowerBoundPtr)
    {
        lowerBounds_.reset(new scalarField(getVars().size(), *lowerBoundPtr));
    }

    // Read upper bounds for the design variables, if present
    if (dict_.found("upperBounds"))
    {
        scalarField upperBounds(dict_.get<scalarField>("upperBounds"));
        if (upperBounds.size() != getVars().size())
        {
            FatalErrorInFunction
                << "Inconsistent dimensions for upperBounds ("
                << upperBounds.size()
                << ") and design variables ("
                << getVars().size() << ")"
                << exit(FatalError);
        }
        upperBounds_.reset(new scalarField(upperBounds));
    }
    else if (dict_.found("upperBound"))
    {
        scalar upperBound(dict_.get<scalar>("upperBound"));
        upperBounds_.reset(new scalarField(getVars().size(), upperBound));
    }
    else if (upperBoundPtr)
    {
        upperBounds_.reset(new scalarField(getVars().size(), *upperBoundPtr));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::designVariables::designVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    scalarField(0),
    mesh_(mesh),
    dict_(dict),
    activeDesignVariables_(0),
    oldDesignVariables_(nullptr),
    maxInitChange_(nullptr),
    lowerBounds_(nullptr),
    upperBounds_(nullptr)
{
    // Read max initial change of design variables if present
    if (dict.found("maxInitChange"))
    {
        maxInitChange_.reset(new scalar(dict_.get<scalar>("maxInitChange")));
    }
}


Foam::designVariables::designVariables
(
    fvMesh& mesh,
    const dictionary& dict,
    const label size
)
:
    scalarField(size, Zero),
    mesh_(mesh),
    dict_(dict),
    activeDesignVariables_(0),
    oldDesignVariables_(nullptr),
    maxInitChange_(nullptr),
    lowerBounds_(nullptr),
    upperBounds_(nullptr)
{
    // Read max initial change of design variables if present
    if (dict.found("maxInitChange"))
    {
        maxInitChange_.reset(new scalar(dict_.get<scalar>("maxInitChange")));
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::designVariables> Foam::designVariables::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    if (!dict.found("type"))
    {
        return autoPtr<designVariables>(nullptr);
    }

    const word modelType(dict.get<word>("type"));

    Info<< "designVariables type : " << modelType << endl;

    auto* ctorPtr = designVariablesConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "designVariables",
            modelType,
            *designVariablesConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<designVariables>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool Foam::designVariables::readDict(const dictionary& dict)
{
    dict_ = dict;

    if (dict.found("maxInitChange"))
    {
        maxInitChange_.reset(new scalar(dict_.get<scalar>("maxInitChange")));
    }

    return true;
}


const Foam::scalarField& Foam::designVariables::getVars() const
{
    return *this;
}


Foam::scalarField& Foam::designVariables::getVars()
{
    return *this;
}


void Foam::designVariables::storeDesignVariables()
{
    if (!oldDesignVariables_)
    {
        oldDesignVariables_.reset(new scalarField(getVars().size(), Zero));
    }

    oldDesignVariables_.ref() = getVars();
}


void Foam::designVariables::resetDesignVariables()
{
    DebugInfo
        << "Reseting design variables" << endl;
    getVars() = (oldDesignVariables_());
}


void Foam::designVariables::postProcessSens
(
    scalarField& objectiveSens,
    PtrList<scalarField>& constraintSens,
    const wordList& adjointSolversNames,
    bool isMaster
)
{
    // Does nothing in base
}


void Foam::designVariables::evolveNumber()
{
    // Does nothing in base
}


void Foam::designVariables::setInitialValues()
{
    // Does nothing in base
}


void Foam::designVariables::addFvOptions
(
    const PtrList<primalSolver>& primalSolver,
    const PtrList<adjointSolverManager>& adjointSolverManagers
)
{
    // Does nothing in base
}


Foam::tmp<Foam::scalarField> Foam::designVariables::constraintValues()
{
    return nullptr;
}


Foam::PtrList<Foam::scalarField> Foam::designVariables::constraintDerivatives()
{
    return PtrList<scalarField>();
}


void Foam::designVariables::writeDesignVars()
{
    // Does nothing in base
}


// ************************************************************************* //
