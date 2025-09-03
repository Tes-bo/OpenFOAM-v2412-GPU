/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "cellBox.H"
#include "mapDistribute.H"

namespace Foam
{
namespace processorLODs
{
    defineTypeNameAndDebug(cellBox, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::processorLODs::cellBox::calcSrcBox
(
    const label srcObji
) const
{
    treeBoundBox bb;

    for (const label facei : srcCells_[srcObji])
    {
        bb.add(srcPoints_, srcFaces_[facei]);
    }

    return bb;
}


Foam::treeBoundBox Foam::processorLODs::cellBox::calcTgtBox
(
    const label tgtObji
) const
{
    treeBoundBox bb;

    for (const label facei : tgtCells_[tgtObji])
    {
        bb.add(tgtPoints_, tgtFaces_[facei]);
    }

    return bb;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::processorLODs::cellBox::cellBox
(
    const cellList& srcCells,
    const faceList& srcFaces,
    const UList<point>& srcPoints,
    const cellList& tgtCells,
    const faceList& tgtFaces,
    const UList<point>& tgtPoints,
    const label maxObjectsPerLeaf,
    const label nObjectsOfType,
    const label nRefineIterMax
)
:
    processorLODs::faceBox
    (
        srcFaces,
        srcPoints,
        tgtFaces,
        tgtPoints,
        maxObjectsPerLeaf,
        nObjectsOfType,
        nRefineIterMax
    ),
    srcCells_(srcCells),
    tgtCells_(tgtCells)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute>
Foam::processorLODs::cellBox::map
(
    const mapDistributeBase::layoutTypes constructLayout
)
{
    return createMap(srcCells_.size(), tgtCells_.size(), constructLayout);
}


// ************************************************************************* //
