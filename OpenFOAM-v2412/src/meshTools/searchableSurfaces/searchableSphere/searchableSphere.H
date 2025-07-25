/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::searchableSphere

Description
    Searching on general spheroid.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        type        | sphere                                | selector |
        origin      | The origin (centre) of the sphere     | yes |
        radius      | The (outside) radius/radiii of sphere | yes |
        centre      | Alternative name for 'origin'         | no  |
    \endtable

Note
    The \c radius can be specified as a single \em scalar (for a sphere)
    or a \em vector of three values (for a general spheroid).

    Longer type name : \c searchableSphere

SourceFiles
    searchableSphere.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSphere_H
#define searchableSphere_H

#include "treeBoundBox.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class searchableSphere Declaration
\*---------------------------------------------------------------------------*/

class searchableSphere
:
    public searchableSurface
{
public:

    // Public Types

        //- The type of shape
        enum shapeType : uint8_t
        {
            SPHERE = 0,     //!< Sphere (all components equal)
            OBLATE = 1,     //!< Oblate (major = mezzo > minor)
            PROLATE = 2,    //!< Prolate (major > mezzo = minor)
            GENERAL = 3     //!< General spheroid
        };


private:

    // Private Types

        //- Component order (largest to smallest)
        struct componentOrder
        {
            uint8_t major;      //!< Component with major radius
            uint8_t mezzo;      //!< Component with intermediate radius
            uint8_t minor;      //!< Component with minor radius
            shapeType shape;    //!< The type of shape
        };


    // Private Data

        //- Centre point of the sphere
        const point origin_;

        //- The outer radii of the spheroid
        const vector radii_;

        //- The canonical (sorted) order and shape
        const struct componentOrder order_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Determine sorted order and classify the shape
        inline static componentOrder getOrdering(const vector& radii);

        //- Shift point relative to origin
        //- and scale relative to spheroid dimensions
        inline point scalePoint(const point& p) const
        {
            return point
            (
                (p.x() - origin_.x()) / radii_.x(),
                (p.y() - origin_.y()) / radii_.y(),
                (p.z() - origin_.z()) / radii_.z()
            );
        }

        //- Undo scalePoint(): unscale point and unshift relative to origin
        inline point unscalePoint(const point& p) const
        {
            return point
            (
                p.x() * radii_.x() + origin_.x(),
                p.y() * radii_.y() + origin_.y(),
                p.z() * radii_.z() + origin_.z()
            );
        }


        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        //- Find nearest point on general spheroid.
        //  With some optimization for special shapes
        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        //- Find intersection with general spheroid
        void findLineAll
        (
            const point& start,
            const point& end,
            pointIndexHit& near,
            pointIndexHit& far
        ) const;


        //- No copy construct
        searchableSphere(const searchableSphere&) = delete;

        //- No copy assignment
        void operator=(const searchableSphere&) = delete;


public:

    //- Runtime type information
    TypeName("searchableSphere");


    // Constructors

        //- Construct a sphere from components
        searchableSphere
        (
            const IOobject& io,
            const point& origin,
            const scalar radius
        );

        //- Construct a spheroid from components
        searchableSphere
        (
            const IOobject& io,
            const point& centre,
            const vector& radii
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableSphere
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSphere() = default;


    // Member Functions

    // Geometric

        //- The centre (origin) of the sphere
        const point& centre() const noexcept
        {
             return origin_;
        }

        //- The radius of the sphere, or major radius of the spheroid
        scalar radius() const noexcept
        {
            return radii_[order_.major];
        }

        //- The radii of the spheroid
        const vector& radii() const noexcept
        {
            return radii_;
        }

        //- The type of shape
        enum shapeType shape() const noexcept
        {
            return order_.shape;
        }

        //- A point on the sphere at given location
        //  theta [-pi,pi], phi [0,pi]
        point surfacePoint(const scalar theta, const scalar phi) const;

        //- Surface normal on the sphere at given location
        //  theta [-pi,pi], phi [0,pi]
        vector surfaceNormal(const scalar theta, const scalar phi) const;



    // Searching

        //- Names of regions
        virtual const wordList& regions() const;

        //- Whether supports volume type (below)
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
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const
        {
            return tmp<pointField>::New(1, origin_);
        }

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const
        {
            return coordinates();
        }

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const;


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

            //- Determine type (inside/outside/mixed) for point.
            virtual void getVolumeType
            (
                const pointField& points,
                List<volumeType>& volType
            ) const;


    // Output

        // Implementation for regIOobject. NotImplemented
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
