/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "maxDurationCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(maxDurationCondition, 0);
    addToRunTimeSelectionTable
    (
        runTimeCondition,
        maxDurationCondition,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::maxDurationCondition::
maxDurationCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    duration_(dict.get<scalar>("duration")),
    startTime_(-1),
    initialised_(false),
    resetOnRestart_(dict.getOrDefault("resetOnRestart", false))
{
    if
    (
        !resetOnRestart_
     && conditionDict().readIfPresent("startTime", startTime_))
    {
        initialised_ = true;
    }
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::maxDurationCondition::apply()
{
    if (!active_)
    {
        return true;
    }

    if (!initialised_)
    {
        startTime_ = obr_.time().value();
        initialised_ = true;
    }

    scalar delta = obr_.time().value() - startTime_;
    delta = obr_.time().timeToUserTime(delta);

    Log << "    " << type() << ": " << name_ << nl
        << "        Completed " << delta << " out of " << duration_ << nl;

    return delta >= duration_;
}


void Foam::functionObjects::runTimeControls::maxDurationCondition::write()
{
    if (initialised_)
    {
        conditionDict().set("startTime", startTime_);
    }
}


void Foam::functionObjects::runTimeControls::maxDurationCondition::reset()
{
    initialised_ = false;
}


// ************************************************************************* //
