/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "clock.H"

#include <string>
#include <sstream>
#include <iomanip>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

static const char *monthNames[] =
{
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

time_t Foam::clock::getTime()
{
    return ::time(reinterpret_cast<time_t*>(0));
}


const struct tm Foam::clock::rawDate()
{
    time_t t = getTime();
    struct tm *curr = ::localtime(&t);
    return *curr;
}


std::string Foam::clock::dateTime()
{
    time_t t = getTime();
    struct tm *curr = ::localtime(&t);

    std::ostringstream os;
    os
        << std::setfill('0')
        << std::setw(4) << curr->tm_year + 1900
        << '-' << std::setw(2) << curr->tm_mon + 1
        << '-' << std::setw(2) << curr->tm_mday
        << 'T'
        << std::setw(2) << curr->tm_hour
        << ':' << std::setw(2) << curr->tm_min
        << ':' << std::setw(2) << curr->tm_sec;

    return os.str();
}


std::string Foam::clock::date()
{
    time_t t = getTime();
    struct tm *curr = ::localtime(&t);

    std::ostringstream os;
    os
        << monthNames[curr->tm_mon]
        << ' ' << std::setw(2) << std::setfill('0') << curr->tm_mday
        << ' ' << std::setw(4) << curr->tm_year + 1900;

    return os.str();
}


std::string Foam::clock::clockTime()
{
    time_t t = getTime();
    struct tm *curr = ::localtime(&t);

    std::ostringstream os;
    os
        << std::setfill('0')
        << std::setw(2) << curr->tm_hour
        << ':' << std::setw(2) << curr->tm_min
        << ':' << std::setw(2) << curr->tm_sec;

    return os.str();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clock::clock()
:
    start_(getTime()),
    last_(start_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double Foam::clock::elapsedClockTime() const
{
    last_ = getTime();
    return ::difftime(last_, start_);
}


double Foam::clock::clockTimeIncrement() const
{
    const auto prev(last_);

    last_ = getTime();
    return ::difftime(last_, prev);
}


// ************************************************************************* //
