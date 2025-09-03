/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "multiDimPolyFunctions.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiDimPolyFunctions, 0);
    defineRunTimeSelectionTable(multiDimPolyFunctions, word);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiDimPolyFunctions> Foam::multiDimPolyFunctions::New
(
    const word& multiDimPolyFunctionsType,
    const labelVector& dirs
)
{
    auto* ctorPtr = wordConstructorTable(multiDimPolyFunctionsType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "multiDimPolyFunction",
            multiDimPolyFunctionsType,
            *wordConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<multiDimPolyFunctions>(ctorPtr(dirs));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiDimPolyFunctions::multiDimPolyFunctions(const labelVector& dirs)
:
    nTerms_(-1),
    geomDir_(dirs),
    geomCorrection_(pos0(dirs.x()), pos0(dirs.y()), pos0(dirs.z())),
    coeffs_(),
    termValues_()
{}


// ************************************************************************* //
