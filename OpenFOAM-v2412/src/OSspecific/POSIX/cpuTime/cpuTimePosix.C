/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "cpuTimePosix.H"
#include <unistd.h>

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

// Clock-ticks per second
static const long clockTicks_(sysconf(_SC_CLK_TCK));


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline double Foam::cpuTimePosix::diff(const value_type& a, const value_type& b)
{
    return
    (
        double((a.tms_utime + a.tms_stime) - (b.tms_utime + b.tms_stime))
      / clockTicks_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuTimePosix::value_type::value_type()
{
    update();
}


Foam::cpuTimePosix::cpuTimePosix()
:
    start_(),
    last_(start_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cpuTimePosix::value_type::update()
{
    ::times(this);
}


void Foam::cpuTimePosix::resetCpuTime()
{
    last_.update();
    start_ = last_;
}


void Foam::cpuTimePosix::resetCpuTimeIncrement() const
{
    last_.update();
}


double Foam::cpuTimePosix::elapsedCpuTime() const
{
    last_.update();
    return diff(last_, start_);
}


double Foam::cpuTimePosix::cpuTimeIncrement() const
{
    const value_type prev(last_);
    last_.update();
    return diff(last_, prev);
}


// ************************************************************************* //
