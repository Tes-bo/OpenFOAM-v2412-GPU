/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "meshState.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::meshState::controlsDictName = "controls";
const Foam::word Foam::meshState::meshDictName = "mesh";
const Foam::word Foam::meshState::solverPerformanceDictName = "solver";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::meshState::getBoolEntry
(
    const dictionary& dict,
    const word& keyword
)
{
    bool value = false;
    return dict.readIfPresent(keyword, value) && value;
}


void Foam::meshState::setBoolEntry
(
    dictionary& dict,
    const word& keyword,
    bool on
)
{
    if (on)
    {
        dict.set(keyword, true);
    }
    else
    {
        dict.remove(keyword);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshState::meshState
(
    const word& name,
    const objectRegistry& obr,
    const dictionary* content
)
:
    IOdictionary
    (
        IOobject
        (
            name,
            obr.time().system(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        content
    ),
    prevTimeIndex_(-1)
{
    (void)subDictOrAdd(controlsDictName);
    (void)subDictOrAdd(meshDictName);
    (void)subDictOrAdd(solverPerformanceDictName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshState::reset(const meshState& ms)
{
    dictionary::operator=(ms);
}


Foam::dictionary& Foam::meshState::controlsDict()
{
    return subDict(controlsDictName);
}


const Foam::dictionary& Foam::meshState::controlsDict() const
{
    return subDict(controlsDictName);
}


Foam::dictionary& Foam::meshState::meshDict()
{
    return subDict(meshDictName);
}


const Foam::dictionary& Foam::meshState::meshDict() const
{
    return subDict(meshDictName);
}


Foam::dictionary& Foam::meshState::solverPerformanceDict()
{
    return subDict(solverPerformanceDictName);
}


const Foam::dictionary& Foam::meshState::solverPerformanceDict() const
{
    return subDict(solverPerformanceDictName);
}


bool Foam::meshState::isFirstIteration() const
{
    return getBoolEntry(controlsDict(), "firstIteration");
}


bool Foam::meshState::isFinalIteration() const
{
    return getBoolEntry(controlsDict(), "finalIteration");
}


void Foam::meshState::setFirstIteration(bool on)
{
    return setBoolEntry(controlsDict(), "firstIteration", on);
}


void Foam::meshState::setFinalIteration(bool on)
{
    return setBoolEntry(controlsDict(), "finalIteration", on);
}


// ************************************************************************* //
