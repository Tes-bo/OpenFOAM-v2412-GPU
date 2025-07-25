/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "cv2DControls.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cv2DControls::cv2DControls
(
    const dictionary& controlDict,
    const boundBox& bb
)
:
    motionControl_(controlDict.subDict("motionControl")),
    conformationControl_(controlDict.subDict("surfaceConformation")),

    minCellSize_(motionControl_.get<scalar>("minCellSize")),
    minCellSize2_(Foam::sqr(minCellSize_)),

    maxQuadAngle_(conformationControl_.get<scalar>("maxQuadAngle")),

    nearWallAlignedDist_
    (
        motionControl_.get<scalar>("nearWallAlignedDist") * minCellSize_
    ),
    nearWallAlignedDist2_(Foam::sqr(nearWallAlignedDist_)),

    insertSurfaceNearestPointPairs_
    (
        conformationControl_.get<Switch>
        (
            "insertSurfaceNearestPointPairs"
        )
    ),
    mirrorPoints_
    (
        conformationControl_.get<Switch>
        (
            "mirrorPoints"
        )
    ),
    insertSurfaceNearPointPairs_
    (
        conformationControl_.get<Switch>
        (
            "insertSurfaceNearPointPairs"
        )
    ),

    objOutput_
    (
        motionControl_.getOrDefault<Switch>("objOutput", false)
    ),

    meshedSurfaceOutput_
    (
        motionControl_.getOrDefault<Switch>("meshedSurfaceOutput", false)
    ),

    randomiseInitialGrid_
    (
        conformationControl_.get<Switch>("randomiseInitialGrid")
    ),
    randomPerturbation_
    (
        conformationControl_.get<scalar>("randomPerturbation")
    ),

    maxBoundaryConformingIter_
    (
        conformationControl_.get<label>("maxBoundaryConformingIter")
    ),

    span_
    (
        max(mag(bb.max().x()), mag(bb.min().x()))
      + max(mag(bb.max().y()), mag(bb.min().y()))
    ),
    span2_(Foam::sqr(span_)),

    minEdgeLen_
    (
        conformationControl_.get<scalar>("minEdgeLenCoeff") * minCellSize_
    ),
    minEdgeLen2_(Foam::sqr(minEdgeLen_)),

    maxNotchLen_
    (
        conformationControl_.get<scalar>("maxNotchLenCoeff") * minCellSize_
    ),
    maxNotchLen2_(Foam::sqr(maxNotchLen_)),

    minNearPointDist_
    (
        conformationControl_.get<scalar>("minNearPointDistCoeff")*minCellSize_
    ),
    minNearPointDist2_(Foam::sqr(minNearPointDist_)),

    ppDist_
    (
        conformationControl_.get<scalar>("pointPairDistanceCoeff")*minCellSize_
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cv2DControls::write(Ostream& os) const
{
    const auto oldLevel = os.indentLevel(1);
    os.precision(2);
    os.flags(ios_base::scientific);

    os << nl << "Outputting CV2D Mesher controls:" << nl
       << token::BEGIN_BLOCK << nl
       << indent << "minCellSize2_         : " << minCellSize2_ << nl
       << indent << "span_ / span2_        : " << span_ << " / " << span2_ << nl
       << indent << "maxNotchLen2_         : " << maxNotchLen2_ << nl
       << indent << "minNearPointDist2_    : " << minNearPointDist2_ << nl
       << indent << "nearWallAlignedDist2_ : " << nearWallAlignedDist2_ << nl
       << indent << "ppDist_               : " << ppDist_ << nl
       << indent << "minEdgeLen2_          : " << minEdgeLen2_ << nl
       << token::END_BLOCK << endl;

    os.indentLevel(oldLevel);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const cv2DControls& s)
{
    s.write(os);
    return os;
}



// ************************************************************************* //
