/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "AMIInterpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::AMIInterpolation> Foam::AMIInterpolation::New
(
    const word& modelName,
    const dictionary& dict,
    const bool reverseTarget
)
{
    DebugInfo << "Selecting model " << modelName << endl;

    auto* ctorPtr = dictConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            typeName,
            modelName,
            *dictConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<AMIInterpolation>(ctorPtr(dict, reverseTarget));
}


Foam::autoPtr<Foam::AMIInterpolation> Foam::AMIInterpolation::New
(
    const word& modelName,
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection
)
{
    DebugInfo << "Selecting model " << modelName << endl;

    auto* ctorPtr = componentConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            typeName,
            modelName,
            *componentConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<AMIInterpolation>
    (
        ctorPtr
        (
            requireMatch,
            reverseTarget,
            lowWeightCorrection
        )
    );
}

// ************************************************************************* //
