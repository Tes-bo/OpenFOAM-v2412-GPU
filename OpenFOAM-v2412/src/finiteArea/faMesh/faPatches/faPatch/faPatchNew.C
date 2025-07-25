/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "faPatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::faPatch> Foam::faPatch::New
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
{
    word patchType(dict.get<word>("type"));
    // Not needed:  dict.readIfPresent("geometricType", patchType);

    return faPatch::New(patchType, name, dict, index, bm);
}


Foam::autoPtr<Foam::faPatch> Foam::faPatch::New
(
    const word& patchType,
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
{
    DebugInFunction << "Constructing faPatch" << endl;

    auto* ctorPtr = dictionaryConstructorTable(patchType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "faPatch",
            patchType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<faPatch>(ctorPtr(name, dict, index, bm, patchType));
}


// ************************************************************************* //
