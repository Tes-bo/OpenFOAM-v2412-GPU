/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

template<class DataType>
inline Foam::PointIntegrateData<DataType>::PointIntegrateData()
:
    valid_(false)
{}


template<class DataType>
inline Foam::PointIntegrateData<DataType>::PointIntegrateData
(
    const DataType& data
)
:
    valid_(true),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DataType>
template<class TrackingData>
inline bool Foam::PointIntegrateData<DataType>::valid(TrackingData& td) const
{
    return valid_;
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointIntegrateData<DataType>::sameGeometry
(
    const PointIntegrateData<DataType>&,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class DataType>
template<class TrackingData>
inline void Foam::PointIntegrateData<DataType>::leaveDomain
(
    const polyPatch& patch,
    const label patchPointi,
    const point& pos,
    TrackingData& td
)
{}


template<class DataType>
template<class TrackingData>
inline void Foam::PointIntegrateData<DataType>::enterDomain
(
    const polyPatch& patch,
    const label patchPointi,
    const point& pos,
    TrackingData& td
)
{}


template<class DataType>
template<class TrackingData>
inline void Foam::PointIntegrateData<DataType>::transform
(
    const tensor& rotTensor,
    TrackingData& td
)
{
    this->data_ = Foam::transform(rotTensor, this->data_);
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointIntegrateData<DataType>::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const label edgeI,
    const PointIntegrateData<DataType>& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Update point from an edge
    if (!valid_)
    {
        if (!edgeInfo.valid_)
        {
            FatalErrorInFunction<< "edgeInfo:" << edgeInfo << exit(FatalError);
        }
        this->operator=(edgeInfo);
        return true;
    }

    return false;
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointIntegrateData<DataType>::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const PointIntegrateData<DataType>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Update point from coupled point
    if (!valid_)
    {
        if (!newPointInfo.valid_)
        {
            FatalErrorInFunction<< "newPointInfo:" << newPointInfo
                << exit(FatalError);
        }
        this->operator=(newPointInfo);
        return true;
    }

    return false;
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointIntegrateData<DataType>::updatePoint
(
    const PointIntegrateData<DataType>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid_)
    {
        if (!newPointInfo.valid_)
        {
            FatalErrorInFunction<< "newPointInfo:" << newPointInfo
                << exit(FatalError);
        }
        this->operator=(newPointInfo);
        return true;
    }

    return false;
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointIntegrateData<DataType>::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointI,
    const PointIntegrateData<DataType>& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid_)
    {
        if (!pointInfo.valid(td))
        {
            FatalErrorInFunction<< "problem: invalid point:"
                << mesh.points()[pointI]
                << " data:" << pointInfo
                << exit(FatalError);
        }
        this->data_ = pointInfo.data_ + td.edgeData_[edgeI];
        this->valid_ = true;
        return true;
    }

    return false;
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointIntegrateData<DataType>::equal
(
    const PointIntegrateData<DataType>& pi,
    TrackingData& td
) const
{
    if (!valid_)
    {
        return false;
    }
    else if (!pi.valid_)
    {
        FatalErrorInFunction << "pi:" << pi
            << exit(FatalError);
        return false;
    }
    else
    {
        return this->data_ == pi.data_;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class DataType>
inline bool Foam::PointIntegrateData<DataType>::operator==
(
    const PointIntegrateData<DataType>& rhs
) const
{
    return this->data_ == rhs.data_;
}


template<class DataType>
inline bool Foam::PointIntegrateData<DataType>::operator!=
(
    const PointIntegrateData<DataType>& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class DataType>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PointIntegrateData<DataType>& rhs
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        return os << rhs.valid_ << token::SPACE << rhs.data();
    }
    else
    {
        return os << rhs.valid_ << rhs.data();
    }
}


template<class DataType>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    PointIntegrateData<DataType>& rhs
)
{
    return is >> rhs.valid_ >> rhs.data_;
}


// ************************************************************************* //
