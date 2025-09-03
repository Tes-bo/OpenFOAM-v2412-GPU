/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021,2023 OpenCFD Ltd.
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

#include "GAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::GAMGInterfaceField> Foam::GAMGInterfaceField::New
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
{
    const word coupleType(fineInterface.interfaceFieldType());

    auto* ctorPtr = lduInterfaceFieldConstructorTable(coupleType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "GAMGInterfaceField",
            coupleType,
            *lduInterfaceFieldConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<GAMGInterfaceField>(ctorPtr(GAMGCp, fineInterface));
}


Foam::autoPtr<Foam::GAMGInterfaceField> Foam::GAMGInterfaceField::New
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
{
    const word coupleType(GAMGCp.type());

    auto* ctorPtr = lduInterfaceConstructorTable(coupleType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "GAMGInterfaceField",
            coupleType,
            *lduInterfaceConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<GAMGInterfaceField>(ctorPtr(GAMGCp, doTransform, rank));
}


Foam::autoPtr<Foam::GAMGInterfaceField> Foam::GAMGInterfaceField::New
(
    const word& patchFieldType,
    const GAMGInterface& GAMGCp,
    Istream& is
)
{
    auto* ctorPtr = IstreamConstructorTable(patchFieldType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "GAMGInterfaceField",
            patchFieldType,
            *IstreamConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<GAMGInterfaceField>(ctorPtr(GAMGCp, is));
}


// ************************************************************************* //
