/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::faceAreaIntersect

Description
    Face intersection class
    - calculates intersection area by sub-dividing face into triangles
      and cutting

SourceFiles
    faceAreaIntersect.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaIntersect_H
#define faceAreaIntersect_H

#include "pointField.H"
#include "FixedList.H"
#include "plane.H"
#include "face.H"
#include "triangle.H"
#include "Enum.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class faceAreaIntersect Declaration
\*---------------------------------------------------------------------------*/

class faceAreaIntersect
{
public:

    enum triangulationMode
    {
        tmFan,
        tmMesh
    };

    static const Enum<triangulationMode> triangulationModeNames_;


private:

    // Private data

        //- Reference to the points of face A
        const pointField& pointsA_;

        //- Reference to the points of face B
        const pointField& pointsB_;

        //- Triangle decomposition of face A
        const DynamicList<face>& trisA_;

        //- Triangle decomposition of face B
        const DynamicList<face>& trisB_;

        //- Flag to reverse B faces
        const bool reverseB_;

        //- Flag to cache the final triangulation
        bool cacheTriangulation_;

        //- Final triangulation
        mutable DynamicList<triPoints> triangles_;


    // Static data members

        //- Tolerance
        static scalar tol;


    // Private Member Functions

        //- Get triPoints from face
        inline triPoints getTriPoints
        (
            const pointField& points,
            const face& f,
            const bool reverse
        ) const;

        //- Set triPoints into tri list
        inline void setTriPoints
        (
            const point& a,
            const point& b,
            const point& c,
            label& count,
            FixedList<triPoints, 10>& tris
        ) const;

        //- Return point of intersection between plane and triangle edge
        inline point planeIntersection
        (
            const FixedList<scalar, 3>& d,
            const triPoints& t,
            const label negI,
            const label posI
        ) const;

        //- Slice triangle with plane and generate new cut sub-triangles
        void triSliceWithPlane
        (
            const triPoints& tri,
            const plane& pln,
            FixedList<triPoints, 10>& tris,
            label& nTris,
            const scalar len
        ) const;

        //- Return area of intersection of triangles src and tgt
        void triangleIntersect
        (
            const triPoints& src,
            const point& tgt0,
            const point& tgt1,
            const point& tgt2,
            const vector& n,
            scalar& area,
            vector& centroid
        ) const;


public:

    // Constructors

        //- Construct from components
        faceAreaIntersect
        (
            const pointField& pointsA,
            const pointField& pointsB,
            const DynamicList<face>& trisA,
            const DynamicList<face>& trisB,
            const bool reverseB = false,
            const bool cacheTriangulation = false
        );


    // Public Member Functions

        //- Fraction of local length scale to use as intersection tolerance
        inline static scalar& tolerance();

        //- Triangulate a face using the given triangulation mode
        static void triangulate
        (
            const face& f,
            const pointField& points,
            const triangulationMode& triMode,
            faceList& faceTris
        );

        //- Const access to the cacheTriangulation flag
        inline bool cacheTriangulation() const;

        //- Const access to the triangulation
        inline const DynamicList<triPoints> triangles() const;

        //- Non-const access to the triangulation
        inline DynamicList<triPoints>& triangles();

        //- Decompose face into triangle fan
        static inline void triangleFan
        (
            const face& f,
            DynamicList<face>& faces
        );

        //- Return area of intersection of faceA with faceB and effective
        //- face centre
        void calc
        (
            const face& faceA,
            const face& faceB,
            const vector& n,
            scalar& area,
            vector& centroid
        ) const;

        //- Return area of intersection of faceA with faceB
        bool overlaps
        (
            const face& faceA,
            const face& faceB,
            const vector& n,
            const scalar threshold
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceAreaIntersectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
