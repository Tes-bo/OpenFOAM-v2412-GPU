/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Class
    Foam::PointData

Description
    Variant of pointEdgePoint with some transported additional data. Templated
    on the transported data type.

SourceFiles
    PointDataI.H
    PointData.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PointData_H
#define Foam_PointData_H

#include "pointEdgePoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class DataType> class PointData;

template<class DataType>
Ostream& operator<<(Ostream&, const PointData<DataType>&);
template<class DataType>
Istream& operator>>(Istream&, PointData<DataType>&);

/*---------------------------------------------------------------------------*\
                           Class PointData Declaration
\*---------------------------------------------------------------------------*/

template<class DataType>
class PointData
:
    public pointEdgePoint
{
    // Private Data

        //- Additional transported data
        DataType data_;


public:

    // Constructors

        //- Default construct
        PointData() = default;

        //- Construct from origin, distance and data
        inline PointData
        (
            const point& origin,
            const scalar distSqr,
            const DataType& data
        );


    // Member Functions

    // Access

        const DataType& data() const noexcept { return data_; }

        DataType& data() noexcept { return data_; }


    // Member Operators

        //- Test for equality
        inline bool operator==(const PointData<DataType>&) const;

        //- Test for inequality
        inline bool operator!=(const PointData<DataType>&) const;


    // IOstream Operators

        friend Ostream& operator<< <DataType>
        (
            Ostream&,
            const PointData<DataType>&
        );
        friend Istream& operator>> <DataType>
        (
            Istream&,
            PointData<DataType>&
        );


    // Wave Methods

        // Needed by MeshWave

            //- Apply rotation matrix to the data
            template<class TrackingData>
            inline void transform
            (
                const tensor& rotTensor,
                TrackingData& td
            );

            //- Influence of edge on point
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const label edgeI,
                const PointData<DataType>& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const PointData<DataType>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const PointData<DataType>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge.
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointI,
                const PointData<DataType>& pointInfo,
                const scalar tol,
                TrackingData& td
            );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if data type is contiguous
template<class DataType>
struct is_contiguous<PointData<DataType>> : is_contiguous<DataType> {};

//- Contiguous scalar only when data type is also scalar
template<class DataType>
struct is_contiguous_scalar<PointData<DataType>>
:
    is_contiguous_scalar<DataType>
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PointDataI.H"

#ifdef NoRepository
#   include "PointData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
