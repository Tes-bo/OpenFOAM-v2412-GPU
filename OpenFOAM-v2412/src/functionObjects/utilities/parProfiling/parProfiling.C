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

#include "parProfiling.H"
#include "profilingPstream.H"
#include "Pstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(parProfiling, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        parProfiling,
        dictionary
    );

} // End namespace functionObject
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::parProfiling
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    reportLevel_(0)
{
    dict.readIfPresent("detail", reportLevel_);
    profilingPstream::enable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::~parProfiling()
{
    profilingPstream::disable();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::parProfiling::report()
{
    if (UPstream::parRun() && UPstream::nProcs() > 1)
    {
        Info<< nl;
        profilingPstream::report(reportLevel_);
    }
}


bool Foam::functionObjects::parProfiling::execute()
{
    report();
    return true;
}


bool Foam::functionObjects::parProfiling::write()
{
    return true;
}


bool Foam::functionObjects::parProfiling::end()
{
    profilingPstream::disable();
    return true;
}


// ************************************************************************* //
