/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subCycleTime::subCycleTime(Time& runTime, const label nCycles)
:
    time_(runTime),
    index_(0),
    total_(nCycles)
{
    // Could avoid 0 or 1 nCycles here on construction
    if (nCycles > 1)
    {
        time_.subCycle(nCycles);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::subCycleTime::~subCycleTime()
{
    endSubCycle();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::subCycleTime::status() const
{
    return (index_ <= total_);
}


bool Foam::subCycleTime::end() const
{
    return (index_ > total_);  // or !(status())
}


void Foam::subCycleTime::endSubCycle()
{
    if (total_ > 1)
    {
        time_.endSubCycle();
    }

    // If called manually, ensure status() will return false

    index_ = total_ + 1;
}


bool Foam::subCycleTime::loop()
{
    const bool active = status();

    if (active)
    {
        operator++();
    }

    return active;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::subCycleTime& Foam::subCycleTime::operator++()
{
    if (total_ > 1)
    {
        time_++;
    }

    index_++;

    // Register index change with Time, in case someone wants this information
    time_.subCycleIndex(index_);

    return *this;
}


Foam::subCycleTime& Foam::subCycleTime::operator++(int)
{
    return operator++();
}


// ************************************************************************* //
