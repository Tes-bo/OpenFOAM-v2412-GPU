/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "instant.H"
#include "Time.H"
#include "Pair.H"
#include "UList.H"
#include <cstdlib>  // std::atof
#include <utility>  // std::move

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::instant::typeName = "instant";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::instant::findStart
(
    const UList<instant>& times,
    const scalar timeVal
)
{
    for (label i = 0; i < times.size(); ++i)
    {
        if (timeVal <= times[i].value())
        {
            return i;
        }
    }
    return 0;
}


Foam::Pair<Foam::label> Foam::instant::findRange
(
    const UList<instant>& times,
    const scalar timeVal,
    const label start
)
{
    Pair<label> range(start, -1);  // lower/upper

    for (label i = start+1; i < times.size(); ++i)
    {
        if (timeVal < times[i].value())
        {
            break;
        }
        else
        {
            range.first() = i;
        }
    }

    if (range.first() < 0 || range.first() >= times.size())
    {
        // Invalid
        return Pair<label>(-1, -1);
    }

    if (range.first() < times.size()-1)
    {
        // Upper part of range within bounds
        range.second() = range.first()+1;
    }

    return range;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::instant::instant(scalar timeValue)
:
    Instant<word>(timeValue, Time::timeName(timeValue))
{}


Foam::instant::instant(const word& timeName)
:
    Instant<word>(0, timeName)
{
    value() = std::atof(name().c_str());
}


Foam::instant::instant(word&& timeName)
:
    Instant<word>(0, std::move(timeName))
{
    value() = std::atof(name().c_str());
}


// ************************************************************************* //
