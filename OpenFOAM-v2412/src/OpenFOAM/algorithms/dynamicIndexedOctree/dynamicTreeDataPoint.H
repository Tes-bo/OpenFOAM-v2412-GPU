/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::dynamicTreeDataPoint

Description
    Holds (reference to) pointField. Encapsulation of data needed for
    octree searches.
    Used for searching for nearest point. No bounding boxes around points.
    Only overlaps and calcNearest are implemented, rest makes little sense.

    Optionally works on subset of points.

SourceFiles
    dynamicTreeDataPoint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dynamicTreeDataPoint_H
#define Foam_dynamicTreeDataPoint_H

#include "pointField.H"
#include "treeBoundBox.H"
#include "line.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class dynamicIndexedOctree;

/*---------------------------------------------------------------------------*\
                           Class dynamicTreeDataPoint Declaration
\*---------------------------------------------------------------------------*/

class dynamicTreeDataPoint
{
    // Private Data

        //- The point field
        const DynamicList<point>& points_;


public:

    // Declare name of the class
    ClassNameNoDebug("dynamicTreeDataPoint");


    // Constructors (non-caching)

        //- Construct from List. Holds reference!
        explicit dynamicTreeDataPoint(const DynamicList<point>& points);


    // Member Functions

        //- Object dimension == 0 (point element)
        int nDim() const noexcept { return 0; }

        //- Return bounding box for the specified point indices
        treeBoundBox bounds(const labelUList& indices) const;


    // Access

        //- The original point field
        const DynamicList<point>& points() const noexcept { return points_; }

        //TDB //- The subset of point ids to use
        ///const labelList& pointLabels() const noexcept { labelList::null(); }

        //- Use a subset of points
        bool useSubset() const noexcept { return false; }

        //- Is the effective point field empty?
        bool empty() const noexcept { return points_.empty(); }

        //- The size of the effective point field
        label size() const noexcept { return points_.size(); }

        //- Map from shape index to original (non-subset) point label
        label objectIndex(label index) const noexcept { return index; }

        //- Point at specified shape index
        const point& operator[](label index) const { return points_[index]; }

        //- Point at specified shape index
        const point& centre(label index) const { return points_[index]; }

        //- Representative point cloud
        const DynamicList<point>& centres() const noexcept { return points_; }


    // Search

        //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
        //  Only makes sense for closed surfaces.
        volumeType getVolumeType
        (
            const dynamicIndexedOctree<dynamicTreeDataPoint>&,
            const point&
        ) const;

        //- Does (bb of) shape at index overlap bb
        bool overlaps
        (
            const label index,
            const treeBoundBox& searchBox
        ) const;

        //- Check if any point on shape is inside sphere.
        bool overlaps
        (
            const label index,
            const point& centre,
            const scalar radiusSqr
        ) const;

        //- Calculates nearest (to sample) point in shape.
        //  Returns actual point and distance (squared)
        void findNearest
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& nearestIndex,
            point& nearestPoint
        ) const;

        //- Calculates nearest (to line) point in shape.
        //  Returns point and distance (squared)
        void findNearest
        (
            const labelUList& indices,
            const linePointRef& ln,

            treeBoundBox& tightest,
            label& minIndex,
            point& linePoint,
            point& nearestPoint
        ) const;

        //- Calculate intersection of shape with ray.
        //  Sets result accordingly
        bool intersects
        (
            const label index,
            const point& start,
            const point& end,
            point& result
        ) const
        {
            NotImplemented;
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
