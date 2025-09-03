/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "extractFromNone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{
    addNamedToRunTimeSelectionTable
    (
        method,
        extractFromNone,
        dictionary,
        none
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::extractFromNone::extractFromNone
(
    const dictionary& dict
)
:
    method()
{
    // A "noneCoeffs" sub-dictionary doesn't make much sense.

    dict.readIfPresent("includedAngle", includedAngle_);
    dict.readIfPresent("geometricTestOnly", geometricTestOnly_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceFeatures>
Foam::surfaceFeaturesExtraction::extractFromNone::features
(
    const triSurface& surf
) const
{
    return autoPtr<surfaceFeatures>::New(surf);
}


// ************************************************************************* //
