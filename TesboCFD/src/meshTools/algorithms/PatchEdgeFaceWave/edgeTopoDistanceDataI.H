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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class PrimitivePatchType>
inline
Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::edgeTopoDistanceData()
:
    distance_(-1),
    data_()
{}


template<class Type, class PrimitivePatchType>
inline
Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::edgeTopoDistanceData
(
    const label distance,
    const Type& data
)
:
    distance_(distance),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::valid
(
    TrackingData& td
) const
{
    return distance_ != -1;
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline void Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::transform
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::updateEdge
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const label edgeI,
    const label facei,
    const edgeTopoDistanceData<Type, PrimitivePatchType>& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From face to edge
    if (distance_ == -1)
    {
        data_ = faceInfo.data_;
        distance_ = faceInfo.distance_ + 1;
        return true;
    }

    return false;
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::updateEdge
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const edgeTopoDistanceData<Type, PrimitivePatchType>& edgeInfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    // From edge to edge (e.g. coupled edges)
    if (distance_ == -1)
    {
        this->operator=(edgeInfo);
        return true;
    }

    return false;
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::updateFace
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const label facei,
    const label edgeI,
    const edgeTopoDistanceData<Type, PrimitivePatchType>& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From edge to face
    if (distance_ == -1)
    {
        this->operator=(edgeInfo);
        return true;
    }

    return false;
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::equal
(
    const edgeTopoDistanceData<Type, PrimitivePatchType>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class PrimitivePatchType>
inline bool Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::operator==
(
    const Foam::edgeTopoDistanceData<Type, PrimitivePatchType>& rhs
) const
{
    return distance() == rhs.distance() && data() == rhs.data();
}


template<class Type, class PrimitivePatchType>
inline bool Foam::edgeTopoDistanceData<Type, PrimitivePatchType>::operator!=
(
    const Foam::edgeTopoDistanceData<Type, PrimitivePatchType>& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
