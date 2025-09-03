/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::patchEdgeFaceRegion::update
(
    const patchEdgeFaceRegion& w2,
    const scalar tol,
    TrackingData& td
)
{
    if (!w2.valid(td))
    {
        FatalErrorInFunction
            << "problem." << abort(FatalError);
    }

    if (w2.region_ == -2 || region_ == -2)
    {
        // Blocked edge/face
        return false;
    }

    if (!valid(td))
    {
        // current not yet set so use any value
        operator=(w2);
        return true;
    }
    else
    {
        if (w2.region_ < region_)
        {
            operator=(w2);
            return true;
        }
        else
        {
            return false;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::patchEdgeFaceRegion::patchEdgeFaceRegion()
:
    region_(-1)
{}


inline Foam::patchEdgeFaceRegion::patchEdgeFaceRegion
(
    const label regioni
)
:
    region_(regioni)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::patchEdgeFaceRegion::valid(TrackingData& td) const
{
    return region_ != -1;
}


template<class TrackingData>
inline void Foam::patchEdgeFaceRegion::transform
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::patchEdgeFaceRegion::updateEdge
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const label edgeI,
    const label facei,
    const patchEdgeFaceRegion& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(faceInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceRegion::updateEdge
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const patchEdgeFaceRegion& edgeInfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    return update(edgeInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceRegion::updateFace
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const label facei,
    const label edgeI,
    const patchEdgeFaceRegion& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(edgeInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceRegion::equal
(
    const patchEdgeFaceRegion& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::patchEdgeFaceRegion::operator==
(
    const patchEdgeFaceRegion& rhs
) const
{
    return region_ == rhs.region_;
}


inline bool Foam::patchEdgeFaceRegion::operator!=
(
    const patchEdgeFaceRegion& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
