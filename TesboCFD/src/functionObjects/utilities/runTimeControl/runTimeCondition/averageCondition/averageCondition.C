/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "averageCondition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(averageCondition, 0);
    addToRunTimeSelectionTable(runTimeCondition, averageCondition, dictionary);
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::averageCondition::averageCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    valueAverageBase(name, obr_, dict, state, false),
    nIterStartUp_(dict.getOrDefault<label>("nIterStartUp", 10)),
    iter_(-1)
{
    dictionary& conditionDict = this->conditionDict();

    readState(conditionDict);

    conditionDict.readIfPresent("iter", iter_);
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::averageCondition::apply()
{
    if (!active_)
    {
        return true;
    }

    bool running = iter_ > nIterStartUp_;

    ++iter_;

    dictionary& conditionDict = this->conditionDict();


    Info<< incrIndent;
    running = valueAverageBase::calculate(conditionDict) && running;
    Info<< decrIndent;

    return running;
}


void Foam::functionObjects::runTimeControls::averageCondition::write()
{
    dictionary& conditionDict = this->conditionDict();

    valueAverageBase::writeState(conditionDict);

    conditionDict.set("iter", iter_);
}


void Foam::functionObjects::runTimeControls::averageCondition::reset()
{
    valueAverageBase::resetState(this->conditionDict());

    iter_ = 0;
}


// ************************************************************************* //
