/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "PtrListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::addPatches(const polyBoundaryMesh& pbm)
{
    // Set boundary patches
    fvPatchList& patches = *this;

    patches.resize_null(pbm.size());

    forAll(patches, patchi)
    {
        patches.set(patchi, fvPatch::New(pbm[patchi], *this));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m
)
:
    fvPatchList(),
    mesh_(m)
{}


Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m,
    const polyBoundaryMesh& pbm
)
:
    fvPatchList(),
    mesh_(m)
{
    addPatches(pbm);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvBoundaryMesh::indices
(
    const wordRe& matcher,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(matcher, useGroups);
}


Foam::labelList Foam::fvBoundaryMesh::indices
(
    const wordRes& matcher,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(matcher, useGroups);
}


Foam::labelList Foam::fvBoundaryMesh::indices
(
    const wordRes& select,
    const wordRes& ignore,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(select, ignore, useGroups);
}


Foam::label Foam::fvBoundaryMesh::findPatchID(const word& patchName) const
{
    if (patchName.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, patchName);
}


void Foam::fvBoundaryMesh::movePoints()
{
    fvPatchList& patches = *this;

    for (fvPatch& p : patches)
    {
        p.initMovePoints();
    }

    for (fvPatch& p : patches)
    {
        p.movePoints();
    }
}


Foam::UPtrList<const Foam::labelUList>
Foam::fvBoundaryMesh::faceCells() const
{
    const fvPatchList& patches = *this;

    UPtrList<const labelUList> list(patches.size());

    forAll(list, patchi)
    {
        list.set(patchi, &patches[patchi].faceCells());
    }

    return list;
}


Foam::lduInterfacePtrsList Foam::fvBoundaryMesh::interfaces() const
{
    const fvPatchList& patches = *this;

    lduInterfacePtrsList list(patches.size());

    forAll(list, patchi)
    {
        const lduInterface* lduPtr = isA<lduInterface>(patches[patchi]);

        if (lduPtr)
        {
            list.set(patchi, lduPtr);
        }
    }

    return list;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::readUpdate(const polyBoundaryMesh& pbm)
{
    addPatches(pbm);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::fvPatch& Foam::fvBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


Foam::fvPatch& Foam::fvBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


// ************************************************************************* //
