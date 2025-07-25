/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "symmetryFaPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symmetryFaPatch, 0);
    addToRunTimeSelectionTable(faPatch, symmetryFaPatch, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::symmetryFaPatch::makeCorrVecs(vectorField& cv) const
{
    // Non-orthogonal correction not allowed
    cv = Zero;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::symmetryFaPatch::symmetryFaPatch
(
    const word& name,
    const labelUList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm,
    const label nbrPolyPatchi,
    const word& patchType
)
:
    faPatch(name, edgeLabels, index, bm, nbrPolyPatchi, patchType)
{}


Foam::symmetryFaPatch::symmetryFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm,
    const word& patchType
)
:
    faPatch(name, dict, index, bm, patchType)
{
    if (ngbPolyPatchIndex() < 0)
    {
        FatalErrorInFunction
            << "Neighbour polyPatch index is not specified for faPatch "
            << this->name() << exit(FatalError);
    }
}


Foam::symmetryFaPatch::symmetryFaPatch
(
    const symmetryFaPatch& p,
    const faBoundaryMesh& bm
)
:
    faPatch(p, bm)
{}


Foam::symmetryFaPatch::symmetryFaPatch
(
    const symmetryFaPatch& p,
    const faBoundaryMesh& bm,
    const label index,
    const labelUList& edgeLabels,
    const label nbrPolyPatchi
)
:
    faPatch(p, bm, index, edgeLabels, nbrPolyPatchi)
{}


// ************************************************************************* //
