/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "adiosTime.H"
#include "adiosReader.H"
#include "TimeState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::adiosFoam::adiosTime::attrType
>
Foam::adiosFoam::adiosTime::attrNames
{
    { attrType::INDEX, "index" },
    { attrType::VALUE, "value" },
    { attrType::DT, "deltaT" },
    { attrType::DT0, "deltaT0" },
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosFoam::adiosTime::adiosTime()
:
    index_(-1),
    value_(VGREAT),
    deltaT_(0),
    deltaT0_(0)
{}


Foam::adiosFoam::adiosTime::adiosTime(const TimeState& t)
:
    index_(t.timeIndex()),
    value_(t.timeOutputValue()),
    deltaT_(t.deltaT().value()),
    deltaT0_(t.deltaT0().value())
{}



Foam::adiosFoam::adiosTime::adiosTime(adiosReader& r)
:
    adiosTime()
{
    if (r.isGood())
    {
        index_ = r.getIntAttribute
        (
            adiosFoam::timeAttribute/attrNames[INDEX]
        );

        value_ = r.getDoubleAttribute
        (
            adiosFoam::timeAttribute/attrNames[VALUE]
        );

        deltaT_ = r.getDoubleAttribute
        (
            adiosFoam::timeAttribute/attrNames[DT]
        );

        deltaT0_ = r.getDoubleAttribute
        (
            adiosFoam::timeAttribute/attrNames[DT0]
        );
    }
}


// ************************************************************************* //
