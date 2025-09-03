/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchDistMethod, 0);
    defineRunTimeSelectionTable(patchDistMethod, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethod::patchDistMethod
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    mesh_(mesh),
    patchIDs_(patchIDs)
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::patchDistMethod>
Foam::patchDistMethod::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const word& defaultPatchDistMethod
)
{
    word modelType(defaultPatchDistMethod);

    // The "method" entry - mandatory if no default was provided
    dict.readEntry
    (
        "method",
        modelType,
        keyType::LITERAL,
        (
            modelType.empty()
          ? IOobjectOption::MUST_READ : IOobjectOption::READ_IF_PRESENT
        )
    );

    Info<< "Selecting patchDistMethod " << modelType << endl;
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "patchDistMethod",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, mesh, patchIDs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchDistMethod::~patchDistMethod()
{}


// ************************************************************************* //
