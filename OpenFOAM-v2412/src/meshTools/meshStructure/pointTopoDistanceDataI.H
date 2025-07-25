/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::pointTopoDistanceData<Type>::pointTopoDistanceData()
:
    distance_(-1),
    data_()
{}


template<class Type>
inline Foam::pointTopoDistanceData<Type>::pointTopoDistanceData
(
    const label distance,
    const Type& data
)
:
    distance_(distance),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<class TrackingData>
inline bool Foam::pointTopoDistanceData<Type>::valid(TrackingData& td) const
{
    return distance_ != -1;
}


// No geometric data so never any problem on cyclics
template<class Type>
template<class TrackingData>
inline bool Foam::pointTopoDistanceData<Type>::sameGeometry
(
    const pointTopoDistanceData<Type>&,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::pointTopoDistanceData<Type>::leaveDomain
(
    const polyPatch& patch,
    const label patchPointi,
    const point& coord,
    TrackingData& td
)
{}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::pointTopoDistanceData<Type>::transform
(
    const tensor& rotTensor,
    TrackingData& td
)
{}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::pointTopoDistanceData<Type>::enterDomain
(
    const polyPatch& patch,
    const label patchPointi,
    const point& coord,
    TrackingData& td
)
{}


// Update this with information from connected edge
template<class Type>
template<class TrackingData>
inline bool Foam::pointTopoDistanceData<Type>::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const pointTopoDistanceData<Type>& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        distance_ = edgeInfo.distance_ + 1;
        data_ = edgeInfo.data_;
        return true;
    }

    return false;
}


// Update this with new information on same point
template<class Type>
template<class TrackingData>
inline bool Foam::pointTopoDistanceData<Type>::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointTopoDistanceData<Type>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        distance_ = newPointInfo.distance_;
        data_ = newPointInfo.data_;
        return true;
    }

    return false;
}


// Update this with new information on same point. No extra information.
template<class Type>
template<class TrackingData>
inline bool Foam::pointTopoDistanceData<Type>::updatePoint
(
    const pointTopoDistanceData<Type>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        distance_ = newPointInfo.distance_;
        data_ = newPointInfo.data_;
        return true;
    }

    return false;
}


// Update this with information from connected point
template<class Type>
template<class TrackingData>
inline bool Foam::pointTopoDistanceData<Type>::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const pointTopoDistanceData<Type>& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        distance_ = pointInfo.distance_;
        data_ = pointInfo.data_;
        return true;
    }

    return false;
}


template<class Type>
template<class TrackingData>
inline bool Foam::pointTopoDistanceData<Type>::equal
(
    const pointTopoDistanceData<Type>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::pointTopoDistanceData<Type>::operator==
(
    const Foam::pointTopoDistanceData<Type>& rhs
) const
{
    return distance() == rhs.distance() && data() == rhs.data();
}


template<class Type>
inline bool Foam::pointTopoDistanceData<Type>::operator!=
(
    const Foam::pointTopoDistanceData<Type>& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
