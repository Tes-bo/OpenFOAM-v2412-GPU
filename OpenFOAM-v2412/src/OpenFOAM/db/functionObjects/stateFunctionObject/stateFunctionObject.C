/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "stateFunctionObject.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionObjects::stateFunctionObject::resultsName_ =
    "results";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::functionObjects::properties&
Foam::functionObjects::stateFunctionObject::stateDict() const
{
    return time_.functionObjects().propsDict();
}


Foam::functionObjects::properties&
Foam::functionObjects::stateFunctionObject::stateDict()
{
    return const_cast<functionObjects::properties&>
    (
        time_.functionObjects().propsDict()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stateFunctionObject::stateFunctionObject
(
    const word& name,
    const Time& runTime
)
:
    timeFunctionObject(name, runTime)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary& Foam::functionObjects::stateFunctionObject::propertyDict()
{
    return stateDict().getObjectDict(name());
}


void Foam::functionObjects::stateFunctionObject::clearTrigger()
{
    return stateDict().clearTrigger();
}


Foam::label Foam::functionObjects::stateFunctionObject::getTrigger() const
{
    return stateDict().getTrigger();
}


bool Foam::functionObjects::stateFunctionObject::setTrigger
(
    const label triggeri
)
{
    return stateDict().setTrigger(triggeri);
}


bool Foam::functionObjects::stateFunctionObject::foundProperty
(
    const word& entryName
) const
{
    return stateDict().foundObjectProperty(name(), entryName);
}


bool Foam::functionObjects::stateFunctionObject::getDict
(
    const word& entryName,
    dictionary& dict
) const
{
    return stateDict().getObjectDict(name(), entryName, dict);
}


bool Foam::functionObjects::stateFunctionObject::getObjectDict
(
    const word& objectName,
    const word& entryName,
    dictionary& dict
) const
{
    return stateDict().getObjectDict(objectName, entryName, dict);
}


bool Foam::functionObjects::stateFunctionObject::getObjectResultDict
(
    const word& objectName,
    dictionary& dict
) const
{
    return stateDict().getObjectResultDict(objectName, dict);
}


Foam::word Foam::functionObjects::stateFunctionObject::resultType
(
    const word& entryName
) const
{
    return stateDict().objectResultType(name(), entryName);
}


Foam::word Foam::functionObjects::stateFunctionObject::objectResultType
(
    const word& objectName,
    const word& entryName
) const
{
    return stateDict().objectResultType(objectName, entryName);
}


Foam::wordList
Foam::functionObjects::stateFunctionObject::objectResultEntries() const
{
    return stateDict().objectResultEntries(name());
}


Foam::wordList
Foam::functionObjects::stateFunctionObject::objectResultEntries
(
    const word& objectName
) const
{
    return stateDict().objectResultEntries(objectName);
}


void Foam::functionObjects::stateFunctionObject::writeResultEntries
(
    Ostream& os
) const
{
    return stateDict().writeResultEntries(name(), os);
}


void Foam::functionObjects::stateFunctionObject::writeResultEntries
(
    const word& objectName,
    Ostream& os
) const
{
    return stateDict().writeResultEntries(objectName, os);
}


void Foam::functionObjects::stateFunctionObject::writeAllResultEntries
(
    Ostream& os
) const
{
    return stateDict().writeAllResultEntries(os);
}


// ************************************************************************* //
