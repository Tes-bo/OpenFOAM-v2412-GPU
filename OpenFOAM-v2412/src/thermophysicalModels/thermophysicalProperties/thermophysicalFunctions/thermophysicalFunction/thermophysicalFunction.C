/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "thermophysicalFunction.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalFunction, 0);
    defineRunTimeSelectionTable(thermophysicalFunction, Istream);
    defineRunTimeSelectionTable(thermophysicalFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalFunction> Foam::thermophysicalFunction::New
(
    Istream& is
)
{
    DebugInFunction << "Constructing thermophysicalFunction" << endl;

    const word functionType(is);

    auto* ctorPtr = IstreamConstructorTable(functionType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "thermophysicalFunction",
            functionType,
            *IstreamConstructorTablePtr_
        ) << abort(FatalError);
    }

    return autoPtr<thermophysicalFunction>(ctorPtr(is));
}


Foam::autoPtr<Foam::thermophysicalFunction> Foam::thermophysicalFunction::New
(
    const dictionary& dict
)
{
    DebugInFunction << "Constructing thermophysicalFunction" << endl;

    const word functionType(dict.get<word>("functionType"));

    auto* ctorPtr = dictionaryConstructorTable(functionType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "thermophysicalFunction",
            functionType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<thermophysicalFunction>(ctorPtr(dict));
}


// ************************************************************************* //
