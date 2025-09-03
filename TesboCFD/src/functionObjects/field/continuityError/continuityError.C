/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "continuityError.H"
#include "volFields.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(continuityError, 0);
    addToRunTimeSelectionTable(functionObject, continuityError, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::continuityError::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Continuity error");

    writeCommented(os, "Time");
    writeCommented(os, "Local");
    writeCommented(os, "Global");
    writeCommented(os, "Cumulative");

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::continuityError::continuityError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    phiName_("phi"),
    cumulative_(getProperty<scalar>("cumulative"))
{
    if (read(dict))
    {
        writeFileHeader(file());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::continuityError::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.readIfPresent("phi", phiName_);

        return true;
    }

    return false;
}


bool Foam::functionObjects::continuityError::execute()
{
    return true;
}


bool Foam::functionObjects::continuityError::write()
{
    const auto* phiPtr = mesh_.cfindObject<surfaceScalarField>(phiName_);

    if (!phiPtr)
    {
        WarningInFunction
            << "Unable to find flux field " << phiName_
            << endl;

        return false;
    }

    const volScalarField error(fvc::div(*phiPtr));
    const scalar deltaT = mesh_.time().deltaTValue();

    const scalar local = deltaT*mag(error)().weightedAverage(mesh_.V()).value();
    const scalar global = deltaT*error.weightedAverage(mesh_.V()).value();
    cumulative_ += global;

    Ostream& os = file();

    writeCurrentTime(os);

    os  << local << tab
        << global << tab
        << cumulative_ << endl;

    Log << type() << " " << name() <<  " write:" << nl
        << "    local = " << local << nl
        << "    global = " << global << nl
        << "    cumulative = " << cumulative_ << nl
        << endl;

    setResult("local", local);
    setResult("global", global);
    setResult("cumulative", cumulative_);

    setProperty<scalar>("cumulative", cumulative_);

    return true;
}


// ************************************************************************* //
