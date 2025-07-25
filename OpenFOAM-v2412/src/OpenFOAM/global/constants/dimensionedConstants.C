/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

Note
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * File Static Data  * * * * * * * * * * * * * //

namespace Foam
{
    dictionary* dimensionedConstantsPtr_(nullptr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::dimensionedConstants()
{
    return debug::switchSet
    (
        "DimensionedConstants",
        dimensionedConstantsPtr_
    );
}


Foam::dimensionedScalar Foam::dimensionedConstant
(
    const word& group,
    const word& varName
)
{
    dictionary& dict = dimensionedConstants();

    // Check that the entries exist.
    // Note: should make FatalError robust instead!

    if (!dict.found("unitSet"))
    {
        std::cerr
            << "Cannot find unitSet in dictionary " << dict.name()
            << std::endl;
    }

    const word unitSetCoeffs(dict.get<word>("unitSet") + "Coeffs");

    const dictionary* unitDictPtr = dict.findDict(unitSetCoeffs);

    if (!unitDictPtr)
    {
        std::cerr
            << "Cannot find " << unitSetCoeffs << " in dictionary "
            << dict.name() << std::endl;
    }

    return dimensionedScalar(varName, unitDictPtr->subDict(group));
}


// ************************************************************************* //
