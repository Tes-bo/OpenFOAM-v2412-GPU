/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::plane

Description
    Geometric class that creates a 3D plane and can return the intersection
    point between a line and the plane.

    Construction from a dictionary is driven by the \c planeType.
    If \c planeType is missing, \c pointAndNormal is used and the
    \c pointAndNormalDict becomes optional.

    For \c planeType as \c pointAndNormal :
    \verbatim
    pointAndNormalDict
    {
        point   <point>;   // basePoint (1612 and earlier)
        normal  <vector>;  // normalVector (1612 and earlier)
    }
    \endverbatim

    For \c planeType as \c embeddedPoints :
    \verbatim
    embeddedPointsDict
    {
        point1  <point>;
        point2  <point>;
        point3  <point>;
    }
    \endverbatim

    For \c planeType with \c planeEquation coefficients
    \f$ ax + by + cz + d = 0 \f$ :

    \verbatim
    planeEquationDict
    {
       a   <scalar>;
       b   <scalar>;
       c   <scalar>;
       d   <scalar>;
    }
    \endverbatim

SourceFiles
    planeI.H
    plane.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_plane_H
#define Foam_plane_H

#include "point.H"
#include "scalarList.H"
#include "line.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                            Class plane Declaration
\*---------------------------------------------------------------------------*/

class plane
{
public:

    // Public Data Types

        //- Side of the plane
        enum side
        {
            FRONT = 1,  //!< The front (positive normal) side of the plane
            BACK  = -1, //!< The back (negative normal) side of the plane
            NORMAL = FRONT,   //!< Alias for FRONT
            FLIP   = BACK     //!< Alias for BACK
        };


    // Public Classes

        //- A reference point and direction
        class ray
        {
            point  point_;
            vector dir_;

        public:

            //- Construct from reference point and direction
            ray(const point& p, const vector& dir)
            :
                point_(p),
                dir_(dir)
            {}

            //- Return the reference point
            const point& refPoint() const noexcept { return point_; }

            //- Return the direction
            const vector& dir() const noexcept { return dir_; }
        };


private:

    // Private Data

        //- The unit normal of the plane
        vector normal_;

        //- The origin or reference point for the plane
        point origin_;


    // Private Member Functions

        //- Normalise normal_ and emit error if its mag is less than VSMALL
        //  Optionally pass as test only, without normalisation.
        void makeUnitNormal
        (
            const char * const caller,
            const bool notTest = true
        );

        //- Calculates point and normal given plane coefficients.
        void calcFromCoeffs
        (
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d,
            const char * const caller
        );

        //- Calculates point and normal vector given three points
        //- Normal vector determined using right hand rule
        void calcFromEmbeddedPoints
        (
            const point& point1,
            const point& point2,
            const point& point3,
            const char * const caller
        );


public:

    // Constructors

        //- Construct zero-initialised.
        inline plane();

        //- Construct from normal vector through the origin.
        //  The vector is normalised to a unit vector on input.
        explicit plane(const vector& normalVector);

        //- Construct from normal vector and point in plane.
        //  By default, the vector is normalised to a unit vector on input.
        plane
        (
            const point& originPoint,
            const vector& normalVector,
            const bool doNormalise = true
        );

        //- Construct from three points in plane
        plane
        (
            const point& point1,
            const point& point2,
            const point& point3
        );

        //- Construct from coefficients for the plane equation:
        //- ax + by + cz + d = 0
        explicit plane(const UList<scalar>& coeffs);

        //- Construct from coefficients for the plane equation:
        //- ax + by + cz + d = 0
        explicit plane(const FixedList<scalar,4>& coeffs);

        //- Construct from dictionary
        explicit plane(const dictionary& dict);

        //- Construct from Istream. Assumes (normal) (origin) input.
        explicit plane(Istream& is);


    // Member Functions

        //- The plane unit normal
        inline const vector& normal() const noexcept;

        //- The plane base point
        inline const point& origin() const noexcept;

        //- The plane base point, for modification
        inline point& origin() noexcept;

        //- Flip the plane by reversing the normal
        inline void flip();

        //- Return coefficients for the plane equation:
        //- ax + by + cz + d = 0
        FixedList<scalar, 4> planeCoeffs() const;

        //- Return nearest point in the plane for the given point
        inline point nearestPoint(const point& p) const;

        //- Return distance (magnitude) from the given point to the plane
        inline scalar distance(const point& p) const;

        //- Return distance from the given point to the plane
        inline scalar signedDistance(const point& p) const;

        //- Return cut coefficient for plane and line defined by
        //- origin and direction
        scalar normalIntersect(const point& pnt0, const vector& dir) const;

        //- Return cut coefficient for plane and ray
        scalar normalIntersect(const ray& r) const
        {
            return normalIntersect(r.refPoint(), r.dir());
        }

        //- Return the cutting point between the plane and
        //- a line passing through the supplied points
        template<class PointType, class PointRef>
        scalar lineIntersect(const line<PointType, PointRef>& l) const
        {
            return normalIntersect(l.start(), l.vec());
        }

        //- Return the cutting line between this plane and another.
        //  Returned as direction vector and point line goes through.
        ray planeIntersect(const plane& plane2) const;

        //- Return the cutting point between this plane and two other planes
        point planePlaneIntersect
        (
            const plane& plane2,
            const plane& plane3
        )
        const;

        //- Return a point somewhere on the plane, a distance from the base
        point somePointInPlane(const scalar dist = 1e-3) const;

        //- Return the side of the plane that the point is on.
        //  If the point is on the plane, then returns FRONT.
        //  \return
        //  - +1 (FRONT): above or on plane (front-side)
        //  - -1 (BACK):  below plane (back-side)
        inline side whichSide(const point& p) const;

        //- The sign for the side of the plane that the point is on.
        //  Uses the supplied tolerance for rounding around zero.
        //  \return
        //  -  0: on plane
        //  - +1 (FRONT): above plane (front-side)
        //  - -1 (BACK):  below plane (back-side)
        inline int sign(const point& p, const scalar tol = SMALL) const;

        //- Mirror the supplied point in the plane. Return the mirrored point.
        point mirror(const point& p) const;

        //- Write to dictionary
        void writeDict(Ostream& os) const;


    // Housekeeping

        //- The plane base point (same as origin)
        const point& refPoint() const noexcept { return origin_; }

        //- Same as whichSide()
        side sideOfPlane(const point& p) const { return whichSide(p); }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Write plane normal, origin
Ostream& operator<<(Ostream& os, const plane& pln);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Test for equality of origin and normal
inline bool operator==(const plane& a, const plane& b);

//- Test for inequality of origin or normal
inline bool operator!=(const plane& a, const plane& b);

//- Compare origin
inline bool operator<(const plane& a, const plane& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "planeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
