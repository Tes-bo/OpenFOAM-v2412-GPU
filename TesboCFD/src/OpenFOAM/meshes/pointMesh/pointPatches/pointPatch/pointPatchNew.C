/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "pointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pointPatch> Foam::pointPatch::New
(
    const word& name,
    const dictionary& dict,
    const label index,
    const pointBoundaryMesh& bm
)
{
    // Similar to polyPatchNew but no support for generic since we want it
    // to fall through to the construct-from-polyPatch
    DebugInFunction << "Constructing pointPatch" << endl;

    const word patchType(dict.lookup("type"));
    //dict.readIfPresent("geometricType", patchType);

    auto* ctorPtr = dictionaryConstructorTable(patchType);

    if (!ctorPtr)
    {
        return nullptr;
    }

    return autoPtr<pointPatch>(ctorPtr(name, dict, index, bm, patchType));
}


// ************************************************************************* //
