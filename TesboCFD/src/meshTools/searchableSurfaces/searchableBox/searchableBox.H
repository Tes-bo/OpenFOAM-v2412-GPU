/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::searchableBox

Description
    Searching on bounding box

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required | Default
        type        | box                               | selector |
        min         | Minimum point for bounding box    | yes   |
        max         | Maximum point for bounding box    | yes   |
        span        | Box width                         | cond'l |
    \endtable

Note
    The order of precedence from the highest to the lowest is
    \c min-max or \c min-span pair

Note
    Longer type name : \c searchableBox

SourceFiles
    searchableBox.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_searchableBox_H
#define Foam_searchableBox_H

#include "searchableSurface.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class searchableBox Declaration
\*---------------------------------------------------------------------------*/

class searchableBox
:
    public searchableSurface,
    public treeBoundBox
{
    // Private Member Data

        mutable wordList regions_;


    // Private Member Functions

        //- Project onto component dir of planePt and update index() (=face)
        void projectOntoCoordPlane
        (
            const direction dir,
            const point& planePt,
            pointIndexHit& info
        ) const;

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        //- Returns miss or hit with face (0..5)
        pointIndexHit findNearest
        (
            const point& bbMid,
            const point& sample,
            const scalar nearestDistSqr
        ) const;


        //- No copy construct
        searchableBox(const searchableBox&) = delete;

        //- No copy assignment
        void operator=(const searchableBox&) = delete;


public:

    //- Runtime type information
    TypeName("searchableBox");


    // Constructors

        //- Construct from components
        searchableBox(const IOobject& io, const treeBoundBox& bb);

        //- Construct from dictionary (used by searchableSurface)
        searchableBox
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableBox() = default;


    // Member Functions

        //- Names of regions
        virtual const wordList& regions() const;

        //- Whether surface supports volume type (below)
        virtual bool hasVolumeType() const
        {
            return true;
        }
        //- What is type of points outside bounds
        virtual volumeType outsideVolumeType() const
        {
            return volumeType::OUTSIDE;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 6;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        // Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const
        {
            return boundBox::overlaps(bb);
        }

        // Single point queries.

            //- Calculate nearest point on surface.
            //  Returns
            //  - bool : any point found nearer than nearestDistSqr
            //  - label: relevant index in surface (=face 0..5)
            //  - point: actual nearest point found
            pointIndexHit findNearest
            (
                const point& sample,
                const scalar nearestDistSqr
            ) const;

            //- Calculate nearest point on edge.
            //  Returns
            //  - bool : any point found nearer than nearestDistSqr
            //  - label: relevant index in surface(=?)
            //  - point: actual nearest point found
            pointIndexHit findNearestOnEdge
            (
                const point& sample,
                const scalar nearestDistSqr
            ) const;

            //- Find nearest to segment.
            //  Returns
            //  - bool : any point found?
            //  - label: relevant index in shapes (=face 0..5)
            //  - point: actual nearest point found
            //  sets:
            //  - tightest  : bounding box
            //  - linePoint : corresponding nearest point on line
            pointIndexHit findNearest
            (
                const linePointRef& ln,
                treeBoundBox& tightest,
                point& linePoint
            ) const;

            //- Find nearest intersection of line between start and end.
            pointIndexHit findLine
            (
                const point& start,
                const point& end
            ) const;

            //- Find any intersection of line between start and end.
            pointIndexHit findLineAny
            (
                const point& start,
                const point& end
            ) const;


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside) for points.
            virtual void getVolumeType
            (
                const pointField& points,
                List<volumeType>& volType
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
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
