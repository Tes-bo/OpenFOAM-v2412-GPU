/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "timeInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        timeInfo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::timeInfo::writeFileHeader(Ostream& os)
{
    writeCommented(os, "Time");
    writeTabbed(os, "cpuTime");
    writeTabbed(os, "clockTime");

    if (perTimeStep_)
    {
        writeTabbed(os, "cpu/step");
        writeTabbed(os, "clock/step");
    }

    os << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeInfo::timeInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    writeFile(time_, name, typeName, dict),
    cpuTime0_(Zero),
    clockTime0_(Zero),
    perTimeStep_(false)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeInfo::read(const dictionary& dict)
{
    timeFunctionObject::read(dict);
    writeFile::read(dict);

    perTimeStep_ = dict.getOrDefault("perTimeStep", false);
    return true;
}


bool Foam::functionObjects::timeInfo::execute()
{
    return true;
}


bool Foam::functionObjects::timeInfo::write()
{
    if (Pstream::master())
    {
        writeCurrentTime(file());

        const scalar cpuTimeNow(time_.elapsedCpuTime());
        const scalar clockTimeNow(time_.elapsedClockTime());

        file()
            << tab << cpuTimeNow
            << tab << clockTimeNow;

        if (perTimeStep_)
        {
            file()
                << tab << (cpuTimeNow - cpuTime0_)
                << tab << (clockTimeNow - clockTime0_);

            cpuTime0_ = cpuTimeNow;
            clockTime0_ = clockTimeNow;
        }

        file() << nl;
    }

    return true;
}


// ************************************************************************* //
