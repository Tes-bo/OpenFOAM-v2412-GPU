/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::PointIntegrateData

Description
    Integrate along selected edges using PointEdgeWave.

SourceFiles
    PointIntegrateDataI.H

\*---------------------------------------------------------------------------*/

#ifndef PointIntegrateData_H
#define PointIntegrateData_H

#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class DataType> class PointIntegrateData;

template<class DataType>
Ostream& operator<<(Ostream&, const PointIntegrateData<DataType>&);
template<class DataType>
Istream& operator>>(Istream&, PointIntegrateData<DataType>&);

/*---------------------------------------------------------------------------*\
                     Class PointIntegrateData declaration
\*---------------------------------------------------------------------------*/

template<class DataType>
class PointIntegrateData
{
    // Private Data

        //- Valid flag
        bool valid_;

        //- Integrated data
        DataType data_;


public:

    //- Class used to pass extra data
    class trackingData
    {
    public:

        UList<DataType>& edgeData_;

        trackingData(UList<DataType>& edgeData)
        :
            edgeData_(edgeData)
        {}
    };


    // Constructors

        //- Default construct
        inline PointIntegrateData();

        //- Construct from data
        inline PointIntegrateData(const DataType& data);


    // Member Functions

        // Access

            const DataType& data() const
            {
                return data_;
            };
            DataType& data()
            {
                return data_;
            };


        // Needed by PointEdgeWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const PointIntegrateData<DataType>&,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Convert origin to relative vector to leaving point
            //- (= point coordinate)
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyPatch& patch,
                const label patchPointi,
                const point& pos,
                TrackingData& td
            );

            //- Convert relative origin to absolute by adding entering point
            template<class TrackingData>
            inline void enterDomain
            (
                const polyPatch& patch,
                const label patchPointi,
                const point& pos,
                TrackingData& td
            );

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
                const PointIntegrateData<DataType>& edgeInfo,
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
                const PointIntegrateData<DataType>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const PointIntegrateData<DataType>& newPointInfo,
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
                const PointIntegrateData<DataType>& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal
            (
                const PointIntegrateData<DataType>&,
                TrackingData& td
            ) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const PointIntegrateData<DataType>&) const;

        //- Test for inequality
        inline bool operator!=(const PointIntegrateData<DataType>&) const;


    // IOstream Operators

        friend Ostream& operator<< <DataType>
        (
            Ostream&,
            const PointIntegrateData<DataType>&
        );
        friend Istream& operator>> <DataType>
        (
            Istream&,
            PointIntegrateData<DataType>&
        );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if the data type is contiguous
template<class DataType>
struct is_contiguous<PointIntegrateData<DataType>>
:
    is_contiguous<DataType>
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PointIntegrateDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
