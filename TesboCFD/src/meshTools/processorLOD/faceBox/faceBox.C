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

#include "faceBox.H"
#include "mapDistribute.H"

namespace Foam
{
namespace processorLODs
{
    defineTypeNameAndDebug(faceBox, 0);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::processorLODs::faceBox::calcSrcBox
(
    const label srcObji
) const
{
    return treeBoundBox(srcPoints_, srcFaces_[srcObji]);  // No reduce
}


Foam::treeBoundBox Foam::processorLODs::faceBox::calcTgtBox
(
    const label tgtObji
) const
{
    return treeBoundBox(tgtPoints_, tgtFaces_[tgtObji]);  // No reduce
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::processorLODs::faceBox::faceBox
(
    const faceList& srcFaces,
    const UList<point>& srcPoints,
    const faceList& tgtFaces,
    const UList<point>& tgtPoints,
    const label maxObjectsPerLeaf,
    const label nObjectsOfType,
    const label nRefineIterMax
)
:
    processorLODs::box(srcPoints, tgtPoints, maxObjectsPerLeaf, nObjectsOfType),
    srcFaces_(srcFaces),
    tgtFaces_(tgtFaces)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute>
Foam::processorLODs::faceBox::map
(
    const mapDistributeBase::layoutTypes constructLayout
)
{
    return createMap(srcFaces_.size(), tgtFaces_.size(), constructLayout);
}


// ************************************************************************* //
