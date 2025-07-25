/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::triangle2D

Description
    2-D triangle and queries

SourceFiles
    triangle2DI.H
    triangle2D.C

\*---------------------------------------------------------------------------*/

#ifndef triangle2D_H
#define triangle2D_H

#include "vector.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class triangle2D Declaration
\*---------------------------------------------------------------------------*/

class triangle2D
:
    public FixedList<vector2D, 3>
{
    // Private Data

        //- Area
        scalar area_;

        //- Helper for intersections
        static FixedList<vector2D, 7> work_;


public:

    static int debug;

    //- Relative tolerance
    static scalar relTol;

    //- Absolute tolerance
    static scalar absTol;

    //- Construct from 3 2-D points
    triangle2D
    (
        const vector2D& a,
        const vector2D& b,
        const vector2D& c,
        const bool orient = false
    );

    //- Construct from 3 3-D points and axes
    triangle2D
    (
        const vector& a3d,
        const vector& b3d,
        const vector& c3d,
        const vector& axis1,
        const vector& axis2,
        const bool orient = false
    );


    // Member Functions

        //- Returns:
        //-   1 if points are ordered in anti-clockwise direction
        //-  -1 if points are ordered in clockwise direction
        //-   0 if the triangle  has collapsed to a line
        inline label order() const;

        //- Write the triangle in OBJ format
        inline static void writeOBJ
        (
            Ostream& os,
            const triangle2D& tri,
            label offset
        );

        //- Return the number of similar points between two triangles
        inline static label nClosePoints
        (
            const triangle2D& triA,
            const triangle2D& triB
        );

        //- Return the signed area
        inline static scalar area
        (
            const vector2D& a,
            const vector2D& b,
            const vector2D& c
        );

        //- Set the intersection between a line and segment
        //- Return true if lines intersect
        inline static bool lineSegmentIntersectionPoint
        (
            const vector2D& lp1,
            const vector2D& lp2,
            const vector2D& sp1,
            const vector2D& sp2,
            vector2D& intersection
        );

        //- Set the intersection between two lines
        //- Return true if lines intersect
        inline static bool lineIntersectionPoint
        (
            const vector2D& a,
            const vector2D& b,
            const vector2D& c,
            const vector2D& d,
            vector2D& intersection
        );

        //- Return true if lines ab and cd intersect
        inline static bool lineIntersects
        (
            const vector2D& a,
            const vector2D& b,
            const vector2D& c,
            const vector2D& d
        );

        //- Snap [this] triangle's points to those of triB if they are within
        //- absTol
        //  Returns the number of snapped points
        label snapClosePoints(const triangle2D& triB);

        //- Return the intersection centre and area
        void interArea
        (
            const triangle2D& triB,
            vector2D& centre,
            scalar& area
        ) const;

        //- Return the intersection area
        scalar interArea(const triangle2D& triB) const;

        //- Return true if triB overlaps
        bool overlaps(const triangle2D& triB) const;

        //- Return the triangle area
        inline scalar area() const noexcept;

        //- Return the triangle centre
        inline vector2D centre() const;

        //- Return true if tri is within this triangle
        inline bool contains(const triangle2D& tri) const;

        //- Return true if triB is the same as this triangle
        inline bool isSame(const triangle2D& triB) const;

        //- Return true if t point p is inside this triangle
        inline bool pointInside(const vector2D& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triangle2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
