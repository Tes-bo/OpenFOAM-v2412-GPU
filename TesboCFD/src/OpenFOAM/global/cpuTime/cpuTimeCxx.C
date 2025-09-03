/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "cpuTimeCxx.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

inline double Foam::cpuTimeCxx::diff(const value_type& a, const value_type& b)
{
    return std::difftime(a.value, b.value) / CLOCKS_PER_SEC;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuTimeCxx::value_type::value_type()
{
    update();
}


Foam::cpuTimeCxx::cpuTimeCxx()
:
    start_(),
    last_(start_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cpuTimeCxx::value_type::update()
{
    value = std::clock();
}


void Foam::cpuTimeCxx::resetCpuTime()
{
    last_.update();
    start_ = last_;
}


void Foam::cpuTimeCxx::resetCpuTimeIncrement() const
{
    last_.update();
}


double Foam::cpuTimeCxx::elapsedCpuTime() const
{
    last_.update();
    return diff(last_, start_);
}


double Foam::cpuTimeCxx::cpuTimeIncrement() const
{
    const value_type prev(last_);
    last_.update();
    return diff(last_, prev);
}


// ************************************************************************* //
