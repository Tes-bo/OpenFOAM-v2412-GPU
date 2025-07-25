/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "engineTime.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::engineTime> Foam::engineTime::New
(
    const word& name,
    const fileName& rootPath,
    const fileName& caseName,
    const fileName& systemName,
    const fileName& constantName,
    const fileName& dictName
)
{
    IFstream engineDictFile("."/constantName/dictName);

    dictionary dict(engineDictFile);

    const word engineType
    (
        dict.getOrDefault<word>("engineType", "crankConRod")
    );

    Info<< "Selecting engine type " << engineType << endl;

    auto* ctorPtr = dictionaryConstructorTable(engineType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "engine",
            engineType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<engineTime>
    (
        ctorPtr
        (
            name,
            rootPath,
            caseName,
            systemName,
            constantName,
            dictName
        )
    );
}


// ************************************************************************* //
