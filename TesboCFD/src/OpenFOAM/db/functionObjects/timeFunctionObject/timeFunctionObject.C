/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "timeFunctionObject.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeFunctionObject::timeFunctionObject
(
    const word& name,
    const Time& runTime
)
:
    functionObject(name),
    time_(runTime)
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::timeFunctionObject::clearOutputObjects
(
    const wordList& objNames
)
{
    objectRegistry& obr = storedObjects();

    for (const word& objName : objNames)
    {
        obr.checkOut(objName);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::objectRegistry&
Foam::functionObjects::timeFunctionObject::storedObjects()
{
    return const_cast<Time&>(time_).functionObjects().storedObjects();
}


const Foam::objectRegistry&
Foam::functionObjects::timeFunctionObject::storedObjects() const
{
    return time_.functionObjects().storedObjects();
}


// ************************************************************************* //
