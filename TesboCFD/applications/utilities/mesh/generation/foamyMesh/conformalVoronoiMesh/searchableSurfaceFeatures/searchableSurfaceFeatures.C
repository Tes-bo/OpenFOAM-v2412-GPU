/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "searchableSurfaceFeatures.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(searchableSurfaceFeatures, 0);
defineRunTimeSelectionTable(searchableSurfaceFeatures, dict);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::searchableSurfaceFeatures>
Foam::searchableSurfaceFeatures::New
(
    const searchableSurface& surface,
    const dictionary& dict
)
{
    const word modelType(surface.type() + "Features");

    auto* ctorPtr = dictConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "searchableSurfaceFeatures",
            modelType,
            *dictConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<searchableSurfaceFeatures>(ctorPtr(surface, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceFeatures::searchableSurfaceFeatures
(
    const searchableSurface& surface,
    const dictionary& dict
)
:
    surface_(surface),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceFeatures::~searchableSurfaceFeatures()
{}


// ************************************************************************* //
