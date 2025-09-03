/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * Global Data * * * * * * * * * * * * * * * //

// Note hard-coded values are more reliable than other alternatives

const Foam::wordList Foam::fieldTypes::basic
({
    "labelField",               //< labelIOField
    "scalarField",              //< scalarIOField
    "vectorField",              //< vectorOField
    "sphericalTensorField",     //< sphericalTensorIOField
    "symmTensorField",          //< symmTensorIOField
    "tensorField"               //< tensorIOField
});


// Commonly used patch field types

const Foam::word Foam::fieldTypes::emptyType
(
    Foam::fieldTypes::emptyTypeName_()
);

const Foam::word Foam::fieldTypes::calculatedType
(
    Foam::fieldTypes::calculatedTypeName_()
);

const Foam::word Foam::fieldTypes::extrapolatedCalculatedType
(
    Foam::fieldTypes::extrapolatedCalculatedTypeName_()
);

const Foam::word Foam::fieldTypes::processorType
(
    Foam::fieldTypes::processorTypeName_()
);

const Foam::word Foam::fieldTypes::zeroGradientType
(
    Foam::fieldTypes::zeroGradientTypeName_()
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
