/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "fakeErrorFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fakeErrorFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        fakeErrorFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fakeErrorFunctionObject::emitError
(
    const char* what
) const
{

    if (ioError_)
    {
        FatalIOError
            << "Error on " << what << " : " << name() << nl
            << exit(FatalIOError);
    }
    else
    {
        FatalError
            << "Error on " << what << " : " << name() << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fakeErrorFunctionObject::fakeErrorFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    ioError_(false),
    constructError_(false),
    executeError_(false),
    writeError_(false)
{
    read(dict);

    if (constructError_)
    {
        emitError("construct");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fakeErrorFunctionObject::read
(
    const dictionary& dict
)
{
    functionObject::read(dict);

    ioError_ = false;
    constructError_ = false;
    executeError_ = false;
    writeError_ = false;

    dict.readIfPresent("ioError", ioError_);
    dict.readIfPresent("constructError", constructError_);
    dict.readIfPresent("executeError", executeError_);
    dict.readIfPresent("writeError", writeError_);

    Log << "Reading : " << name() << nl
        << "    error on construct " << constructError_ << nl
        << "    error on execute() " << executeError_ << nl
        << "    error on write()   " << writeError_ << nl
        << "    using ioerror = " << ioError_ << nl;

    return true;
}


bool Foam::functionObjects::fakeErrorFunctionObject::execute()
{
    if (executeError_)
    {
        emitError("execute()");
    }

    return true;
}


bool Foam::functionObjects::fakeErrorFunctionObject::write()
{
    if (writeError_)
    {
        emitError("write()");
    }

    return true;
}


// ************************************************************************* //
