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

#include "emptyFaPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// Patch name
defineTypeNameAndDebug(emptyFaPatch, 0);

// Add the patch constructor functions to the hash tables
addToRunTimeSelectionTable(faPatch, emptyFaPatch, dictionary);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emptyFaPatch::emptyFaPatch
(
    const word& name,
    const label index,
    const faBoundaryMesh& bm,
    const label nbrPolyPatchi,
    const word& patchType
)
:
    emptyFaPatch(name, labelList(), index, bm, nbrPolyPatchi, patchType)
{}


Foam::emptyFaPatch::emptyFaPatch
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


Foam::emptyFaPatch::emptyFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm,
    const word& patchType
)
:
    faPatch(name, dict, index, bm, patchType)
{}


Foam::emptyFaPatch::emptyFaPatch
(
    const emptyFaPatch& p,
    const faBoundaryMesh& bm
)
:
    faPatch(p, bm)
{}


Foam::emptyFaPatch::emptyFaPatch
(
    const emptyFaPatch& p,
    const faBoundaryMesh& bm,
    const label index,
    const labelUList& edgeLabels,
    const label nbrPolyPatchi
)
:
    faPatch(p, bm, index, edgeLabels, nbrPolyPatchi)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Over-riding the face normals return from the underlying patch
// This is the only piece of info used out of the underlying primitivePatch
// I choose to store it there because it is used in primitive patch operations
// and it should not be duplicated as before.  However, to ensure everything
// in the empty patch is sized to zero, we shall here return a reference to
// a zero-sized field (it does not matter what the field is
//
// const vectorField& emptyFaPatch::edgeNormals() const
// {
//     return faceAreas();
// }


// ************************************************************************* //
