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

#include "LESdelta.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LESdelta, 0);
    defineRunTimeSelectionTable(LESdelta, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESdelta::LESdelta
(
    const word& name,
    const turbulenceModel& turbulence
)
:
    turbulenceModel_(turbulence),
    delta_
    (
        IOobject
        (
            name,
            turbulence.mesh().time().timeName(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence.mesh(),
        dimensionedScalar(word::null, dimLength, SMALL)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LESdelta> Foam::LESdelta::New
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict,
    const word& lookupName
)
{
    const word deltaType(dict.get<word>(lookupName));

    Info<< "Selecting LES " << lookupName << " type " << deltaType << endl;

    auto* ctorPtr = dictionaryConstructorTable(deltaType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "LESdelta",
            deltaType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<LESdelta>(ctorPtr(name, turbulence, dict));
}


Foam::autoPtr<Foam::LESdelta> Foam::LESdelta::New
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict,
    const dictionaryConstructorTableType& additionalConstructors,
    const word& lookupName
)
{
    const word deltaType(dict.get<word>(lookupName));

    Info<< "Selecting LES " << lookupName << " type " << deltaType << endl;

    // First any additional ones
    {
        auto ctorIter = additionalConstructors.cfind(deltaType);

        if (ctorIter.good())
        {
            return autoPtr<LESdelta>(ctorIter.val()(name, turbulence, dict));
        }
    }

    auto* ctorPtr = dictionaryConstructorTable(deltaType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "LESdelta",
            deltaType,
            *dictionaryConstructorTablePtr_
        );

        if (!additionalConstructors.empty())
        {
            FatalIOError
                << " and " << additionalConstructors.sortedToc() << nl;
        }

        FatalIOError
            << exit(FatalIOError);
    }

    return autoPtr<LESdelta>(ctorPtr(name, turbulence, dict));
}


// ************************************************************************* //
