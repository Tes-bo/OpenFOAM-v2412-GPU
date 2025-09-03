/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "polyPatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& patchType,
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
{
    DebugInFunction << "Constructing polyPatch" << endl;

    auto* ctorPtr = wordConstructorTable(patchType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "polyPatch",
            patchType,
            *wordConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<polyPatch>
    (
        ctorPtr
        (
            name,
            size,
            start,
            index,
            bm,
            patchType
        )
    );
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
{
    DebugInFunction << "Constructing polyPatch" << endl;

    word patchType(dict.get<word>("type"));
    dict.readIfPresent("geometricType", patchType);

    return polyPatch::New(patchType, name, dict, index, bm);
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& patchType,
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
{
    DebugInFunction << "Constructing polyPatch" << endl;

    auto* ctorPtr = dictionaryConstructorTable(patchType);

    if (!ctorPtr)
    {
        if (!disallowGenericPolyPatch)
        {
            ctorPtr = dictionaryConstructorTable("genericPatch");
        }

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "polyPatch",
                patchType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }
    }

    return autoPtr<polyPatch>(ctorPtr(name, dict, index, bm, patchType));
}


// ************************************************************************* //
