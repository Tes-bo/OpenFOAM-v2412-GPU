/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "coordinateRotation.H"
#include "axesRotation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(coordinateRotation);
    defineRunTimeSelectionTable(coordinateRotation, dictionary);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::vector Foam::coordinateRotation::findOrthogonal(const vector& axis)
{
    direction maxCmpt = 0;
    scalar maxVal = mag(axis[maxCmpt]);

    for (direction cmpt = 1; cmpt < vector::nComponents; ++cmpt)
    {
        const scalar val = mag(axis[cmpt]);

        if (maxVal < val)
        {
            maxVal  = val;
            maxCmpt = cmpt;
        }
    }

    direction cmpt = ((maxCmpt == vector::nComponents-1) ? 0 : (maxCmpt+1));

    vector dirn(Zero);
    dirn.component(cmpt) = ((axis[maxCmpt] < 0) ? -1 : 1);

    return dirn;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateRotation>
Foam::coordinateRotation::New
(
    const word& modelType,
    const dictionary& dict
)
{
    // Direct dispatch
    // - treat missing modelType as 'axes' (eg, e1/e3 specification)

    if (modelType.empty())
    {
        return autoPtr<coordinateRotation>(new coordinateRotations::axes(dict));
    }

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "rotation",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<coordinateRotation>(ctorPtr(dict));
}


Foam::autoPtr<Foam::coordinateRotation>
Foam::coordinateRotation::New
(
    const dictionary& dict
)
{
    return coordinateRotation::New(dict.get<word>("type"), dict);
}


// ************************************************************************* //
