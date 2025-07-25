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

Class
    Foam::pointTopoDistanceData

Description
    For use with PointEdgeWave. Determines topological distance to
    starting points. Templated on passive transported data.

SourceFiles
    pointTopoDistanceDataI.H
    pointTopoDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef pointTopoDistanceData_H
#define pointTopoDistanceData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
template<class Type> class pointTopoDistanceData;
template<class Type>
Istream& operator>>(Istream&, pointTopoDistanceData<Type>&);
template<class Type>
Ostream& operator<<(Ostream&, const pointTopoDistanceData<Type>&);

/*---------------------------------------------------------------------------*\
                    Class pointTopoDistanceData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointTopoDistanceData
{

    // Protected data

        //- Distance
        label distance_;

        //- Starting data
        Type data_;


public:

    typedef Type dataType;


    // Constructors

        //- Construct null with invalid (-1) for distance, null constructor
        //  for data
        inline pointTopoDistanceData();

        //- Construct from components
        inline pointTopoDistanceData
        (
            const label distance,
            const Type& data
        );


    // Member Functions

        // Access

            inline label distance() const
            {
                return distance_;
            }

            inline const Type& data() const
            {
                return data_;
            }


        // Needed by PointEdgeWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const pointTopoDistanceData<Type>&,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Convert origin to relative vector to leaving point
            //  (= point coordinate)
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

            //- Apply rotation matrix to origin
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
                const label pointi,
                const label edgeI,
                const pointTopoDistanceData<Type>& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const pointTopoDistanceData<Type>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const pointTopoDistanceData<Type>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge.
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointi,
                const pointTopoDistanceData<Type>& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const pointTopoDistanceData<Type>&, TrackingData&)
            const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const pointTopoDistanceData<Type>&) const;
        inline bool operator!=(const pointTopoDistanceData<Type>&) const;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const pointTopoDistanceData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            pointTopoDistanceData<Type>&
        );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if data type is contiguous
template<class Type>
struct is_contiguous<pointTopoDistanceData<Type>> : is_contiguous<Type> {};

//- Data are contiguous label if data type is label
template<class Type>
struct is_contiguous_label<pointTopoDistanceData<Type>> :
    is_contiguous_label<Type> {};

//- Data are contiguous scalar if data type is scalar
template<class Type>
struct is_contiguous_scalar<pointTopoDistanceData<Type>> :
    is_contiguous_scalar<Type>{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointTopoDistanceData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointTopoDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
