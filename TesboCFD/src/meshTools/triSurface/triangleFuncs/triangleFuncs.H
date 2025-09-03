/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::triangleFuncs

Description
    Contains various triangle static functions.

SourceFiles
    triangleFuncs.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_triangleFuncs_H
#define Foam_triangleFuncs_H

#include "pointField.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class treeBoundBox;

/*---------------------------------------------------------------------------*\
                           Class triangleFuncs Declaration
\*---------------------------------------------------------------------------*/

class triangleFuncs
{
    // Private Member Functions

        //- Sets pt to be anywhere on the edge between oppositeSidePt
        //- and thisSidePt depending on both signs.
        // \note Helper function for intersect()
        static void setIntersection
        (
            const point& oppositeSidePt,
            const scalar oppositeSign,
            const point& thisSidePt,
            const scalar thisSign,
            const scalar tol,
            point& pt
        );

public:

    // Static Member Functions

    //- Intersect triangle with parallel edges aligned with axis i0.
    //  Returns true (and intersection in pInter) if any of them intersects
    //  triangle. Used in intersectBb.
    static bool intersectAxesBundle
    (
        const point& V0,
        const point& V10,
        const point& V20,
        const label i0,
        const pointField& origin,
        const scalar maxLength,
        point& pInter
    );

    //- Intersect triangle with bounding box.
    //  \return true if any bounding box faces intersect the triangle,
    //      returns false if the triangle is inside the bounding box
    static bool intersectBb
    (
        const triPointRef& tri,
        const treeBoundBox& cubeBb
    );

    //- Intersect triangle with bounding box.
    //  \return true if any bounding box faces intersect the triangle,
    //      returns false if the triangle is inside the bounding box
    static bool intersectBb
    (
        const point& p0,
        const point& p1,
        const point& p2,
        const treeBoundBox& cubeBb
    );

    //- Intersect triangle with plane.
    //  \return bool and set intersection segment.
    static bool intersect
    (
        const point& va0,
        const point& va10,
        const point& va20,

        const point& basePoint,
        const vector& normal,

        point& pInter0,
        point& pInter1
    );

    //- Intersection of two triangles intersect.
    //  \return bool and set intersection segment.
    static bool intersect
    (
        const point& va0,
        const point& va10,
        const point& va20,

        const point& vb0,
        const point& vb10,
        const point& vb20,

        point& pInter0,
        point& pInter1
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
