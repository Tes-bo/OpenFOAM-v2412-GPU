/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "patchInteractionDataList.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::patchInteractionDataList::patchInteractionDataList()
:
    List<patchInteractionData>(),
    patchGroupIDs_()
{}


Foam::patchInteractionDataList::patchInteractionDataList
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    List<patchInteractionData>(dict.lookup("patches")),
    patchGroupIDs_(this->size())
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const List<patchInteractionData>& items = *this;
    forAllReverse(items, i)
    {
        const wordRe& patchName = items[i].patchName();
        labelList ids = bMesh.indices(patchName);

        if (ids.empty())
        {
            WarningInFunction
                << "Cannot find any patch names matching "
                << patchName << endl;
        }

        patchGroupIDs_[i].transfer(ids);
    }

    // Check that all patches are specified
    DynamicList<word> badPatches;
    for (const polyPatch& pp : bMesh)
    {
        if
        (
            !pp.coupled()
         && !isA<emptyPolyPatch>(pp)
         && applyToPatch(pp.index()) < 0
        )
        {
            badPatches.append(pp.name());
        }
    }

    if (!badPatches.empty())
    {
        FatalErrorInFunction
            << "All patches must be specified when employing local patch "
            << "interaction. Please specify data for patches:" << nl
            << badPatches << nl
            << exit(FatalError);
    }
}


Foam::patchInteractionDataList::patchInteractionDataList
(
    const patchInteractionDataList& pidl
)
:
    List<patchInteractionData>(pidl),
    patchGroupIDs_(pidl.patchGroupIDs_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::patchInteractionDataList::applyToPatch(const label id) const
{
    forAll(patchGroupIDs_, groupi)
    {
        if (patchGroupIDs_[groupi].found(id))
        {
            return groupi;
        }
    }

    return -1;
}


// ************************************************************************* //
