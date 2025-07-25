/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "cellSizeAndAlignmentControl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellSizeAndAlignmentControl, 0);
defineRunTimeSelectionTable(cellSizeAndAlignmentControl, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSizeAndAlignmentControl::cellSizeAndAlignmentControl
(
    const Time& runTime,
    const word& name,
    const dictionary& dict,
    const conformationSurfaces& geometryToConformTo,
    const scalar& defaultCellSize
)
:
    runTime_(runTime),
    defaultCellSize_(defaultCellSize),
    forceInitialPointInsertion_
    (
        dict.getOrDefault<Switch>
        (
            "forceInitialPointInsertion",
            Switch::OFF
        )
    ),
    name_(name)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cellSizeAndAlignmentControl>
Foam::cellSizeAndAlignmentControl::New
(
    const Time& runTime,
    const word& name,
    const dictionary& dict,
    const conformationSurfaces& geometryToConformTo,
    const scalar& defaultCellSize
)
{
    const word modelType(dict.get<word>("type"));

    Info<< indent << "Selecting cellSizeAndAlignmentControl "
        << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "cellSizeAndAlignmentControl",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<cellSizeAndAlignmentControl>
    (
        ctorPtr
        (
            runTime,
            name,
            dict,
            geometryToConformTo,
            defaultCellSize
        )
    );
}


// ************************************************************************* //
