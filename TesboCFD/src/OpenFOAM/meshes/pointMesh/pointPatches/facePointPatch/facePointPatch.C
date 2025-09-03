/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "facePointPatch.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(facePointPatch, 0);
    defineRunTimeSelectionTable(facePointPatch, polyPatch);

    addToRunTimeSelectionTable
    (
        facePointPatch,
        facePointPatch,
        polyPatch
    );

} // End namespace Foam


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::facePointPatch::initGeometry(PstreamBuffers&)
{}


void Foam::facePointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::facePointPatch::initMovePoints(PstreamBuffers&, const pointField&)
{}


void Foam::facePointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::facePointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initGeometry(pBufs);
}


void Foam::facePointPatch::updateMesh(PstreamBuffers&)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::facePointPatch::facePointPatch
(
    const polyPatch& p,
    const pointBoundaryMesh& bm
)
:
    pointPatch(p.name(), p.index(), bm, p.physicalType(), p.inGroups()),
    polyPatch_(p)
{}


Foam::facePointPatch::facePointPatch
(
    const facePointPatch& pp,
    const pointBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const labelUList& reversePointMap
)
:
    pointPatch
    (
        pp.name(),
        index,
        bm,
        pp.physicalType(),
        pp.inGroups()
    ),
    polyPatch_(pp.patch())
{}


// ************************************************************************* //
