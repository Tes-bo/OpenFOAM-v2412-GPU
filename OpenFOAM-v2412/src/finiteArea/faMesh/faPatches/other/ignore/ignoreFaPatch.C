/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "ignoreFaPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// Patch name
defineTypeNameAndDebug(ignoreFaPatch, 0);

// Add the patch constructor functions to the hash tables
addToRunTimeSelectionTable(faPatch, ignoreFaPatch, dictionary);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ignoreFaPatch::ignoreFaPatch
(
    const word& name,
    const label index,
    const faBoundaryMesh& bm,
    const label nbrPolyPatchi,
    const word& patchType
)
:
    ignoreFaPatch(name, labelList(), index, bm, nbrPolyPatchi, patchType)
{}


Foam::ignoreFaPatch::ignoreFaPatch
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


Foam::ignoreFaPatch::ignoreFaPatch
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


Foam::ignoreFaPatch::ignoreFaPatch
(
    const ignoreFaPatch& p,
    const faBoundaryMesh& bm
)
:
    faPatch(p, bm)
{}


Foam::ignoreFaPatch::ignoreFaPatch
(
    const ignoreFaPatch& p,
    const faBoundaryMesh& bm,
    const label index,
    const labelUList& edgeLabels,
    const label nbrPolyPatchi
)
:
    faPatch(p, bm, index, edgeLabels, nbrPolyPatchi)
{}


// ************************************************************************* //
