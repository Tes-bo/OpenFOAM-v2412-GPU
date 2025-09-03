/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "tabulatedWallFunction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::tabulatedWallFunctions::tabulatedWallFunction>
Foam::tabulatedWallFunctions::tabulatedWallFunction::New
(
    const dictionary& dict,
    const polyMesh& mesh
)
{
    const word functionName(dict.get<word>("tabulatedWallFunction"));

    Info<< "Selecting tabulatedWallFunction " << functionName << endl;

    auto* ctorPtr = dictionaryConstructorTable(functionName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "tabulatedWallFunction",
            functionName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<tabulatedWallFunction>(ctorPtr(dict, mesh));
}


// ************************************************************************* //
