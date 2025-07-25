/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
inline Foam::topoDistanceData<Type>::topoDistanceData()
:
    distance_(-1),
    data_()
{}


template<class Type>
inline Foam::topoDistanceData<Type>::topoDistanceData
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
inline bool Foam::topoDistanceData<Type>::valid(TrackingData& td) const
{
    return distance_ != -1;
}


// No geometric data so never any problem on cyclics
template<class Type>
template<class TrackingData>
inline bool Foam::topoDistanceData<Type>::sameGeometry
(
    const polyMesh&,
    const topoDistanceData<Type>&,
    const scalar,
    TrackingData&
) const
{
    return true;
}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::topoDistanceData<Type>::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData&
)
{}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::topoDistanceData<Type>::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData&
)
{}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::topoDistanceData<Type>::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData&
)
{}


// Update cell with neighbouring face information
template<class Type>
template<class TrackingData>
inline bool Foam::topoDistanceData<Type>::updateCell
(
    const polyMesh&,
    const label thisCelli,
    const label neighbourFacei,
    const topoDistanceData<Type>& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    if (distance_ == -1)
    {
        this->operator=(neighbourInfo);
        return true;
    }

    return false;
}


// Update face with neighbouring cell information
template<class Type>
template<class TrackingData>
inline bool Foam::topoDistanceData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const topoDistanceData<Type>& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From cell to its faces.

    if (distance_ == -1)
    {
        distance_ = neighbourInfo.distance_ + 1;
        data_ = neighbourInfo.data_;
        return true;
    }

    return false;
}


// Update face with coupled face information
template<class Type>
template<class TrackingData>
inline bool Foam::topoDistanceData<Type>::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const topoDistanceData<Type>& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From face to face (e.g. coupled faces)
    if (distance_ == -1)
    {
        this->operator=(neighbourInfo);
        return true;
    }

    return false;
}


template<class Type>
template<class TrackingData>
inline bool Foam::topoDistanceData<Type>::equal
(
    const topoDistanceData<Type>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::topoDistanceData<Type>::operator==
(
    const Foam::topoDistanceData<Type>& rhs
) const
{
    return distance() == rhs.distance() && data() == rhs.data();
}


template<class Type>
inline bool Foam::topoDistanceData<Type>::operator!=
(
    const Foam::topoDistanceData<Type>& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
