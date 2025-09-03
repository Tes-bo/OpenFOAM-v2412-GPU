/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "function1Base.H"
#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::objectRegistry* Foam::function1Base::whichDb
(
    const bool useTime
) const noexcept
{
    if (obrPtr_ && useTime)
    {
        return &(obrPtr_->time());
    }

    return obrPtr_;
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::function1Base::function1Base
(
    const word& entryName,
    const objectRegistry* obrPtr
)
:
    refCount(),
    name_(entryName),
    obrPtr_(obrPtr)
{}


Foam::function1Base::function1Base
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    refCount(),
    name_(entryName),
    obrPtr_(obrPtr)
{}


Foam::function1Base::function1Base(const function1Base& rhs)
:
    refCount(),
    name_(rhs.name_),
    obrPtr_(rhs.obrPtr_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// NOTE : do not delete obrPtr_ (no ownership)
Foam::function1Base::~function1Base()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::function1Base::obr() const
{
    if (!obrPtr_)
    {
        FatalErrorInFunction
            << "Object registry not set"
            << abort(FatalError);
    }
    return *obrPtr_;
}


const Foam::Time& Foam::function1Base::time() const
{
    if (!obrPtr_)
    {
        FatalErrorInFunction
            << "Object registry not set"
            << abort(FatalError);
    }

    return obrPtr_->time();
}


bool Foam::function1Base::isTime() const noexcept
{
    return (obrPtr_ && obrPtr_->isTimeDb());
}


void Foam::function1Base::resetDb(const objectRegistry* obrPtr) noexcept
{
    obrPtr_ = obrPtr;
}


void Foam::function1Base::resetDb(const objectRegistry& db) noexcept
{
    obrPtr_ = &db;
}


void Foam::function1Base::userTimeToTime(const Time& t)
{}


// ************************************************************************* //
