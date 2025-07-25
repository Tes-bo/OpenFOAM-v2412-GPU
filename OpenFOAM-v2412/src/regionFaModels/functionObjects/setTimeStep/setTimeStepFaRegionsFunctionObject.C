/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "setTimeStepFaRegionsFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(setTimeStepFaRegionsFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        setTimeStepFaRegionsFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::
setTimeStepFaRegionsFunctionObject::
setTimeStepFaRegionsFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::setTimeStepFaRegionsFunctionObject::adjustTimeStep()
{
    // Wanted timestep
    scalar newDeltaT = regionDeltaT();

    static label index = -1;

    if ((time_.timeIndex() != index) && (newDeltaT < time_.deltaTValue()))
    {
        // Store current time so we don't get infinite recursion (since
        // setDeltaT calls adjustTimeStep() again)
        index = time_.timeIndex();

        // Set time, allow deltaT to be adjusted for writeInterval purposes
        const_cast<Time&>(time_).setDeltaT(newDeltaT, false);

        return true;
    }

    return false;
}


bool Foam::functionObjects::setTimeStepFaRegionsFunctionObject::read
(
    const dictionary& dict
)
{
    if (timeFunctionObject::read(dict))
    {
        // Ensure that adjustTimeStep is active
        if (!time_.controlDict().getOrDefault<bool>("adjustTimeStep", false))
        {
            FatalIOErrorInFunction(dict)
                << "Need to set 'adjustTimeStep' true to allow timestep control"
                << nl
                << exit(FatalIOError);
        }

        return true;
    }

    return false;
}


Foam::scalar Foam::functionObjects::setTimeStepFaRegionsFunctionObject::
regionDeltaT() const
{
    bool adjust = false;
    scalar Co = 0;

    for (const regionFaModel& reg : time_.cobjects<regionFaModel>())
    {
        const scalar regionCo = reg.CourantNumber();
        if (Co < regionCo)
        {
            Co = regionCo;
            adjust = true;
        }
    }

    if (adjust)
    {
        const scalar regionFaMaxCo =
            time_.controlDict().get<scalar>("regionFaMaxCo");

        const scalar maxDeltaTFact = regionFaMaxCo/(Co + SMALL);
        const scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        return deltaTFact*time_.deltaTValue();
    }

    return time_.deltaTValue();
}


bool Foam::functionObjects::setTimeStepFaRegionsFunctionObject::execute()
{
    return true;
}


bool Foam::functionObjects::setTimeStepFaRegionsFunctionObject::write()
{
    return true;
}


// ************************************************************************* //
