/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "initialPointsMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(initialPointsMethod, 0);
defineRunTimeSelectionTable(initialPointsMethod, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::initialPointsMethod::initialPointsMethod
(
    const word& type,
    const dictionary& initialPointsDict,
    const Time& runTime,
    Random& rndGen,
    const conformationSurfaces& geometryToConformTo,
    const cellShapeControl& cellShapeControls,
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
:
    dictionary(initialPointsDict),
    runTime_(runTime),
    rndGen_(rndGen),
    geometryToConformTo_(geometryToConformTo),
    cellShapeControls_(cellShapeControls),
    decomposition_(decomposition),
    detailsDict_(optionalSubDict(type + "Coeffs")),
    minimumSurfaceDistanceCoeffSqr_
    (
        sqr
        (
            initialPointsDict.get<scalar>("minimumSurfaceDistanceCoeff")
        )
    ),
    fixInitialPoints_(initialPointsDict.get<bool>("fixInitialPoints"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::initialPointsMethod> Foam::initialPointsMethod::New
(
    const dictionary& dict,
    const Time& runTime,
    Random& rndGen,
    const conformationSurfaces& geometryToConformTo,
    const cellShapeControl& cellShapeControls,
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
{
    const word modelType(dict.get<word>("initialPointsMethod"));

    Info<< nl << "Selecting initialPointsMethod " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "initialPointsMethod",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return
        autoPtr<initialPointsMethod>
        (
            ctorPtr
            (
                dict,
                runTime,
                rndGen,
                geometryToConformTo,
                cellShapeControls,
                decomposition
            )
        );
}


// ************************************************************************* //
