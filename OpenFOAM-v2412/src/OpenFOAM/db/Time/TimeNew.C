/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "Time.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::Time> Foam::Time::New()
{
    return autoPtr<Time>::New
    (
        fileName("."),  // root-path
        fileName("."),  // case-name
        false,          // No enableFunctionObjects
        false           // No enableLibs
    );
}

// FUTURE?
// Foam::autoPtr<Foam::Time> Foam::Time::New(const Time& runTime)
// {
//     fileName caseDir(runTime.path());
//     caseDir.toAbsolute();
//
//     return autoPtr<Time>::New
//     (
//         caseDir.path(), // root-path
//         caseDir.name(), // case-name
//         false,          // No enableFunctionObjects
//         false           // No enableLibs
//     );
// }


Foam::autoPtr<Foam::Time> Foam::Time::New(const fileName& caseDir)
{
    return autoPtr<Time>::New
    (
        caseDir.path(), // root-path
        caseDir.name(), // case-name
        false,          // No enableFunctionObjects
        false           // No enableLibs
    );
}


Foam::autoPtr<Foam::Time> Foam::Time::New(const argList& args)
{
    return autoPtr<Time>::New
    (
        Time::controlDictName,
        args,
        false,          // No enableFunctionObjects
        false,          // No enableLibs
        IOobjectOption::MUST_READ  // No re-reading
    );
}


Foam::autoPtr<Foam::Time> Foam::Time::NewGlobalTime()
{
    fileName caseDir(argList::envGlobalPath());
    caseDir.toAbsolute();

    return autoPtr<Time>::New
    (
        caseDir.path(), // root-path
        caseDir.name(), // case-name
        false,          // No enableFunctionObjects
        false           // No enableLibs
    );
}


Foam::autoPtr<Foam::Time> Foam::Time::NewGlobalTime(const Time& runTime)
{
    fileName caseDir(runTime.globalPath());
    caseDir.toAbsolute();

    return autoPtr<Time>::New
    (
        caseDir.path(), // root-path
        caseDir.name(), // case-name
        false,          // No enableFunctionObjects
        false           // No enableLibs
    );
}


// ************************************************************************* //
