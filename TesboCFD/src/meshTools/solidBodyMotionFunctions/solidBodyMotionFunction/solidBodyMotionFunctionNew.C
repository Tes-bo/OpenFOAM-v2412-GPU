/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidBodyMotionFunction>
Foam::solidBodyMotionFunction::New
(
    const word& motionType,
    const dictionary& dict,
    const Time& runTime
)
{
    if (motionType.empty())
    {
        return nullptr;
    }

    auto* ctorPtr = dictionaryConstructorTable(motionType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "solidBodyMotionFunction",
            motionType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    Info<< "Selecting solid-body motion function " << motionType << endl;

    return autoPtr<solidBodyMotionFunction>(ctorPtr(dict, runTime));
}


Foam::autoPtr<Foam::solidBodyMotionFunction>
Foam::solidBodyMotionFunction::New
(
    const dictionary& dict,
    const Time& runTime
)
{
    return New
    (
        dict.get<word>("solidBodyMotionFunction", keyType::LITERAL),
        dict,
        runTime
    );
}


Foam::autoPtr<Foam::solidBodyMotionFunction>
Foam::solidBodyMotionFunction::NewIfPresent
(
    const dictionary& dict,
    const Time& runTime
)
{
    word motionType;
    dict.readIfPresent("solidBodyMotionFunction", motionType, keyType::LITERAL);

    if (motionType.empty())
    {
        return nullptr;
    }

    return New(motionType, dict, runTime);
}


// ************************************************************************* //
