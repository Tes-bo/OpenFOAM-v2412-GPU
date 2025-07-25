/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "lumpedPointController.H"
#include "dictionary.H"
#include "labelField.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointController::lumpedPointController() noexcept
:
    pointLabels_()
{}


Foam::lumpedPointController::lumpedPointController
(
    const labelUList& pointLabels
)
:
    pointLabels_(pointLabels)
{}


Foam::lumpedPointController::lumpedPointController
(
    labelList&& pointLabels
)
:
    pointLabels_(std::move(pointLabels))
{}


Foam::lumpedPointController::lumpedPointController
(
    const dictionary& dict
)
:
    pointLabels_(dict.get<labelList>("pointLabels"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointController::remapPointLabels
(
    const label nPoints,
    const Map<label>& originalIds
)
{
    if (originalIds.size())
    {
        for (label& pointi : pointLabels_)
        {
            pointi = originalIds[pointi];
        }
    }

    if (min(pointLabels_) < 0 || max(pointLabels_) >= nPoints)
    {
        FatalErrorInFunction
            << "Point id out-of-range: " << flatOutput(pointLabels_) << nl
            << exit(FatalError);
    }
}


// ************************************************************************* //
