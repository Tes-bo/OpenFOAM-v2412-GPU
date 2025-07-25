/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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

#include "regularisationRadius.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regularisationRadius, 0);
    defineRunTimeSelectionTable(regularisationRadius, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regularisationRadius::regularisationRadius
(
    const fvMesh& mesh,
    const dictionary& dict,
    bool adjustWallThickness
)
:
    mesh_(mesh),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::regularisationRadius> Foam::regularisationRadius::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    bool adjustWallThickness
)
{
    const word modelType = dict.getOrDefault<word>("type", "isotropic");

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    Info<< "regularisationRadius type " << modelType << endl;

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "regularisationRadius",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<regularisationRadius>
    (
        ctorPtr(mesh, dict, adjustWallThickness)
    );
}


// ************************************************************************* //
