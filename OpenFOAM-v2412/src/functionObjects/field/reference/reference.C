/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "reference.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reference, 0);
    addToRunTimeSelectionTable(functionObject, reference, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::reference::calc()
{
    Log << type() << " " << name() << " output:" << nl;

    bool processed = calcType<scalar>();
    processed = processed || calcType<vector>();
    processed = processed || calcType<sphericalTensor>();
    processed = processed || calcType<symmTensor>();
    processed = processed || calcType<tensor>();

    Log << endl;

    return returnReduceOr(processed);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reference::reference
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    localDict_(dict),
    scale_(1)
{
    read(dict);

    // Forcing default result name to include the name of the field
    setResultName(typeName, word::null);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reference::read(const dictionary& dict)
{
    if (fieldExpression::read(dict))
    {
        localDict_ = dict;

        Log << type() << " " << name() << nl
            << "    field: " << fieldName_ << nl;

        if (dict.readIfPresent("scale", scale_))
        {
            Log << "    scale: " << scale_ << nl;
        }

        Log << endl;

        return true;
    }

    return false;
}


// ************************************************************************* //
