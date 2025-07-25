/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "runTimeCondition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(runTimeCondition, 0);
    defineRunTimeSelectionTable(runTimeCondition, dictionary);
}
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::dictionary&
Foam::functionObjects::runTimeControls::runTimeCondition::setConditionDict()
{
    dictionary& propertyDict = state_.propertyDict();

    if (!propertyDict.found(name_))
    {
        propertyDict.add(name_, dictionary());
    }

    return propertyDict.subDict(name_);
}


const Foam::dictionary&
Foam::functionObjects::runTimeControls::runTimeCondition::conditionDict() const
{
    return conditionDict_;
}


Foam::dictionary&
Foam::functionObjects::runTimeControls::runTimeCondition::conditionDict()
{
    return conditionDict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::runTimeCondition::runTimeCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    name_(name),
    obr_(obr),
    state_(state),
    active_(dict.getOrDefault("active", true)),
    conditionDict_(setConditionDict()),
    groupID_(dict.getOrDefault("groupID", -1)),
    log(dict.getOrDefault("log", true))
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

const Foam::word&
Foam::functionObjects::runTimeControls::runTimeCondition::name() const
{
    return name_;
}


bool Foam::functionObjects::runTimeControls::runTimeCondition::active() const
{
    return active_;
}


Foam::label
Foam::functionObjects::runTimeControls::runTimeCondition::groupID() const
{
    return groupID_;
}


// ************************************************************************* //
