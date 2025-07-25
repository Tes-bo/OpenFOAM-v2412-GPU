/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "thermophysicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalProperties, 0);
    defineRunTimeSelectionTable(thermophysicalProperties,);
    defineRunTimeSelectionTable(thermophysicalProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalProperties::thermophysicalProperties(scalar W)
:
    W_(W)
{}


Foam::thermophysicalProperties::thermophysicalProperties(const dictionary& dict)
:
    W_(dict.get<scalar>("W"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const word& name
)
{
    DebugInFunction << "Constructing thermophysicalProperties" << endl;

    auto* ctorPtr = ConstructorTable(name);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "thermophysicalProperties",
            name,
            *ConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<thermophysicalProperties>(ctorPtr());
}


Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const dictionary& dict
)
{
    DebugInFunction << "Constructing thermophysicalProperties" << endl;

    const word& modelType = dict.dictName();

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "thermophysicalProperties",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<thermophysicalProperties>(ctorPtr(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermophysicalProperties::readIfPresent(const dictionary &dict)
{
    dict.readIfPresent("W", W_);
}


void Foam::thermophysicalProperties::writeData(Ostream& os) const
{
    os  << W_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const thermophysicalProperties& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
