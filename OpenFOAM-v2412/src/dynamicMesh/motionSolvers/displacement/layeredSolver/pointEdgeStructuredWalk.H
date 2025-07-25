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

Class
    Foam::pointEdgeStructuredWalk

Description
    Determines length of string of edges walked to point.

SourceFiles
    pointEdgeStructuredWalkI.H
    pointEdgeStructuredWalk.C

\*---------------------------------------------------------------------------*/

#ifndef pointEdgeStructuredWalk_H
#define pointEdgeStructuredWalk_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
class pointEdgeStructuredWalk;

Istream& operator>>(Istream&, pointEdgeStructuredWalk&);
Ostream& operator<<(Ostream&, const pointEdgeStructuredWalk&);

/*---------------------------------------------------------------------------*\
                   Class pointEdgeStructuredWalk Declaration
\*---------------------------------------------------------------------------*/

class pointEdgeStructuredWalk
{
    // Private Data

        //- Starting location
        point point0_;

        //- Previous point
        point previousPoint_;

        //- Sum of distance
        scalar dist_;

        //- Passive data
        vector data_;

        //- Index of passive data (optional)
        label index_;


    // Private Member Functions

        //- Evaluate distance to point.
        template<class TrackingData>
        inline bool update
        (
            const pointEdgeStructuredWalk& w2,
            const scalar tol,
            TrackingData& td
        );

public:

    // Constructors

        //- Default construct
        inline pointEdgeStructuredWalk();

        //- Construct from components
        inline pointEdgeStructuredWalk
        (
            const point& point0,
            const point& previousPoint,
            const scalar dist,
            const vector& data,
            const label index = -1
        );


    // Member Functions

    // Access

        //- The distance information
        scalar dist() const
        {
            return dist_;
        }

        //- Tracking data
        const vector& data() const
        {
            return data_;
        }

        //- Index (if any) associated with data
        label index() const
        {
            return index_;
        }

        //- True if starting point is valid (ie, not point::max)
        inline bool inZone() const;


        // Needed by MeshWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const pointEdgeStructuredWalk&,
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
                const pointEdgeStructuredWalk& edgeInfo,
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
                const pointEdgeStructuredWalk& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const pointEdgeStructuredWalk& newPointInfo,
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
                const pointEdgeStructuredWalk& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal
            (
                const pointEdgeStructuredWalk&,
                TrackingData&
            ) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const pointEdgeStructuredWalk&) const;

        //- Test for inequality
        inline bool operator!=(const pointEdgeStructuredWalk&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pointEdgeStructuredWalk&);
        friend Istream& operator>>(Istream&, pointEdgeStructuredWalk&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for pointEdgeStructuredWalk
template<> struct is_contiguous<pointEdgeStructuredWalk> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointEdgeStructuredWalkI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
