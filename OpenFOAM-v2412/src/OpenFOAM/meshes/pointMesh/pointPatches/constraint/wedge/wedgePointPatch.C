/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "wedgePointPatch.H"
#include "pointConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wedgePointPatch, 0);

    // Add the patch constructor functions to the hash tables
    addToRunTimeSelectionTable
    (
        facePointPatch,
        wedgePointPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wedgePointPatch::wedgePointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    facePointPatch(patch, bm),
    wedgePolyPatch_(refCast<const wedgePolyPatch>(patch))
{}


Foam::wedgePointPatch::wedgePointPatch
(
    const wedgePointPatch& patch,
    const pointBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const labelUList& reversePointMap
)
:
    facePointPatch(patch, bm, index, mapAddressing, reversePointMap),
    wedgePolyPatch_(refCast<const wedgePolyPatch>(patch.patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wedgePointPatch::applyConstraint
(
    const label pointi,
    pointConstraint& pc
) const
{
    pc.applyConstraint(wedgePolyPatch_.n());
}


// ************************************************************************* //
