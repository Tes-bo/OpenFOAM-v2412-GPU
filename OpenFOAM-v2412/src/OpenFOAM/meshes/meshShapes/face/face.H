/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::face

Description
    A face is a list of labels corresponding to mesh vertices.

See also
    Foam::triFace

SourceFiles
    faceI.H
    face.C
    faceIntersection.C
    faceContactSphere.C
    faceAreaInContact.C
    faceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_face_H
#define Foam_face_H

#include "pointField.H"
#include "labelList.H"
#include "edgeList.H"
#include "vectorField.H"
#include "faceListFwd.H"
#include "intersection.H"
#include "pointHit.H"
#include "FixedList.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class face;
class triFace;
template<class T, int SizeMin> class DynamicList;

/*---------------------------------------------------------------------------*\
                           Class face Declaration
\*---------------------------------------------------------------------------*/

class face
:
    public labelList
{
    // Private Member Functions

        //- Construct list of edge vectors for face
        tmp<vectorField> calcEdgeVectors
        (
            const UList<point>& points
        ) const;

        //- Find index of largest internal angle on face
        label mostConcaveAngle
        (
            const UList<point>& points,
            const vectorField& edges,
            scalar& edgeCos
        ) const;

        //- Enumeration listing the modes for split()
        enum splitMode
        {
            COUNTTRIANGLE,  //!< count if split into triangles
            COUNTQUAD,      //!< count if split into triangles and quads
            SPLITTRIANGLE,  //!< split into triangles
            SPLITQUAD       //!< split into  triangles and quads
        };

        //- Split face into triangles or triangles and quads.
        //  Stores results quadFaces[quadI], triFaces[triI]
        //  \return number of new faces created
        label split
        (
            const splitMode mode,
            const UList<point>& points,
            label& triI,
            label& quadI,
            faceList& triFaces,
            faceList& quadFaces
        ) const;


public:

    // Public Typedefs

        //- The proximity classifications
        enum proxType
        {
            NONE = 0,   //!< Unknown proximity
            POINT,      //!< Close to point
            EDGE        //!< Close to edge
        };


    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Default construct
        constexpr face() noexcept = default;

        //- Construct given size, with invalid vertex labels (-1)
        inline explicit face(const label sz);

        //- Copy construct from list of vertex labels
        inline explicit face(const labelUList& list);

        //- Move construct from list of vertex labels
        inline explicit face(labelList&& list);

        //- Copy construct from an initializer list of vertex labels
        inline explicit face(std::initializer_list<label> list);

        //- Copy construct from list of vertex labels
        template<unsigned N>
        inline explicit face(const FixedList<label, N>& list);

        //- Copy construct from subset of vertex labels
        inline face(const labelUList& list, const labelUList& indices);

        //- Copy construct from subset of vertex labels
        template<unsigned N>
        inline face(const labelUList& list, const FixedList<label, N>& indices);

        //- Copy construct from triFace
        face(const triFace& f);

        //- Construct from Istream
        inline face(Istream& is);


    // Member Functions

        //- Collapse face by removing duplicate vertex labels
        //  \return the collapsed size
        label collapse();

        //- Flip the face in-place.
        //  The starting vertex of the original and flipped face are identical.
        void flip();

        //- Return the points corresponding to this face
        inline pointField points(const UList<point>& pts) const;

        //- Centre point of face
        point centre(const UList<point>& points) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type average
        (
            const UList<point>& meshPoints,
            const Field<Type>& fld
        ) const;

        //- The area normal - with magnitude equal to area of face
        vector areaNormal(const UList<point>& p) const;

        //- The unit normal
        inline vector unitNormal(const UList<point>& p) const;

        //- Legacy name for areaNormal()
        //  \deprecated(2018-06) Deprecated for new use
        FOAM_DEPRECATED_FOR(2018-12, "areaNormal() or unitNormal()")
        vector normal(const UList<point>& p) const
        {
            return areaNormal(p); // Legacy definition
        }

        //- Magnitude of face area
        inline scalar mag(const UList<point>& p) const;

        //- The enclosing (bounding) box for the face
        inline Pair<point> box(const UList<point>& pts) const;

        //- Return face with reverse direction
        //  The starting vertex of the original and reverse face are identical.
        face reverseFace() const;


        // Navigation through face vertices

        //- Find local vertex on face for the vertex label, same as find()
        //  \return position in face (0,1,2,...) or -1 if not found.
        inline label which(const label vertex) const;

        //- The vertex on face - identical to operator[], but with naming
        //- similar to nextLabel(), prevLabel()
        inline label thisLabel(const label i) const;

        //- Next vertex on face
        inline label nextLabel(const label i) const;

        //- Previous vertex on face
        inline label prevLabel(const label i) const;


        //- Return the volume swept out by the face when its points move
        scalar sweptVol
        (
            const UList<point>& oldPoints,
            const UList<point>& newPoints
        ) const;

        //- Return the inertia tensor, with optional reference
        //- point and density specification
        tensor inertia
        (
            const UList<point>& p,
            const point& refPt = vector::zero,
            scalar density = 1.0
        ) const;

        //- Return potential intersection with face with a ray starting
        //- at p, direction n (does not need to be normalized)
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. Face-centre is calculated from point average.
        //  For a hit, the distance is signed.  Positive number
        //  represents the point in front of triangle
        //  In case of miss the point is the nearest point on the face
        //  and the distance is the distance between the intersection point
        //  and the original point.
        //  The half-ray or full-ray intersection and the contact
        //  sphere adjustment of the projection vector is set by the
        //  intersection parameters
        pointHit ray
        (
            const point& p,
            const vector& n,
            const UList<point>& meshPoints,
            const intersection::algorithm alg = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        ) const;

        //- Fast intersection with a ray.
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. See triangle::intersection for details.
        pointHit intersection
        (
            const point& p,
            const vector& q,
            const point& ctr,
            const UList<point>& meshPoints,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        //- Return nearest point to face
        pointHit nearestPoint
        (
            const point& p,
            const UList<point>& meshPoints
        ) const;

        //- Return nearest point to face and classify it:
        //  + near point (nearType=POINT, nearLabel=0, 1, 2)
        //  + near edge (nearType=EDGE, nearLabel=0, 1, 2)
        //    Note: edges are counted from starting vertex so
        //    e.g. edge n is from f[n] to f[0], where the face has n + 1
        //    points
        pointHit nearestPointClassify
        (
            const point& p,
            const UList<point>& meshPoints,
            label& nearType,
            label& nearLabel
        ) const;

        //- The sign for the side of the face plane the point is on,
        //- using three evenly distributed face points for the estimated normal.
        //  Uses the supplied tolerance for rounding around zero.
        //  \return
        //  -  0: on plane
        //  - +1: front-side
        //  - -1: back-side
        int sign
        (
            const point& p,
            const UList<point>& points,
            const scalar tol = SMALL
        ) const;

        //- Return contact sphere diameter
        scalar contactSphereDiameter
        (
            const point& p,
            const vector& n,
            const UList<point>& meshPoints
        ) const;

        //- Return area in contact, given the displacement in vertices
        scalar areaInContact
        (
            const UList<point>& meshPoints,
            const scalarField& v
        ) const;

        //- Return number of edges
        inline label nEdges() const noexcept;

        //- Return i-th face edge (forward walk order).
        //  The edge 0 is from [0] to [1]
        //  and edge 1 is from [1] to [2]
        inline Foam::edge edge(const label edgei) const;

        //- Return vector of i-th face edge (forward walk order).
        //  The edge 0 is from [0] to [1]
        //  and edge 1 is from [1] to [2]
        inline vector edge(const label edgei, const UList<point>& pts) const;

        //- Return i-th face edge in reverse walk order.
        //  The rcEdge 0 is from [0] to [n-1]
        //  and rcEdge 1 is from [n-1] to [n-2]
        inline Foam::edge rcEdge(const label edgei) const;

        //- Return vector of i-th face edge in reverse walk order.
        //  The rcEdge 0 is from [0] to [n-1]
        //  and rcEdge 1 is from [n-1] to [n-2]
        inline vector rcEdge(const label edgei, const UList<point>& pts) const;

        //- Return list of edges in forward walk order.
        //  The edge 0 is from [0] to [1]
        //  and edge 1 is from [1] to [2]
        edgeList edges() const;

        //- Return list of edges in reverse walk order.
        //  The rcEdge 0 is from [0] to [n-1]
        //  and rcEdge 1 is from [n-1] to [n-2]
        edgeList rcEdges() const;


    // Searching

        //- Regular contains(pointi) tests
        using labelList::contains;

        //- True if face contains(edge)
        inline bool contains(const Foam::edge& e) const;

        //- Regular find pointi within face
        using labelList::find;

        //- Find the edge within the face.
        //  \return the edge index or -1 if not found
        label find(const Foam::edge& e) const;

        //- Test the edge direction on the face
        //  \return
        //  -  0: edge not found on the face
        //  - +1: forward (counter-clockwise) on the face
        //  - -1: reverse (clockwise) on the face
        int edgeDirection(const Foam::edge& e) const;

        //- Find the longest edge on a face.
        label longestEdge(const UList<point>& pts) const;


        // Face splitting utilities

            //- Number of triangles after splitting
            inline label nTriangles() const noexcept;

            //- Number of triangles after splitting
            label nTriangles(const UList<point>& unused) const;

            //- Split into triangles using existing points.
            //  Result in triFaces[triI..triI+nTri].
            //  Splits intelligently to maximize triangle quality.
            //  \Return number of faces created.
            label triangles
            (
                const UList<point>& points,
                label& triI,
                faceList& triFaces
            ) const;

            //- Split into triangles using existing points.
            //  Append to DynamicList.
            //  \Return number of faces created.
            template<int SizeMin>
            label triangles
            (
                const UList<point>& points,
                DynamicList<face, SizeMin>& triFaces
            ) const;

            //- Number of triangles and quads after splitting
            //  Returns the sum of both
            label nTrianglesQuads
            (
                const UList<point>& points,
                label& nTris,
                label& nQuads
            ) const;

            //- Split into triangles and quads.
            //  Results in triFaces (starting at triI) and quadFaces
            //  (starting at quadI).
            //  Returns number of new faces created.
            label trianglesQuads
            (
                const UList<point>& points,
                label& triI,
                label& quadI,
                faceList& triFaces,
                faceList& quadFaces
            ) const;


        //- True if the face has at least one vertex in common with other
        inline bool connected(const labelUList& other) const;

        //- True if the face has at least one vertex in common with other
        template<unsigned N>
        inline bool connected(const FixedList<label, N>& other) const;

        //- Compare faces
        //  \return
        //  -  0: different
        //  - +1: identical
        //  - -1: same face, but different orientation
        static int compare(const face& a, const face& b);

        //- True if the faces have all the same vertices
        static bool sameVertices(const face& a, const face& b);


    // Member Operators

        //- Increment (offset) vertices by given amount
        inline void operator+=(const label vertexOffset);


    // Hashing

        //- The symmetric hash value for face.
        //  Starts with lowest vertex value and walks in the direction
        //  of the next lowest value.
        static unsigned symmhash_code(const UList<label>& f, unsigned seed=0);

        //- The hash value for face.
        //  Currently hashes as sequential contiguous data,
        //  but could/should be commutative
        inline unsigned hash_code(unsigned seed=0) const
        {
            return Foam::Hasher(this->cdata(), this->size_bytes(), seed);
        }

        //- The symmetric hash value for face.
        //  Starts with lowest vertex value and walks in the direction
        //  of the next lowest value.
        inline unsigned symmhash_code(unsigned seed=0) const
        {
            return face::symmhash_code(*this, seed);
        }

        //- Hashing functor for face
        struct hasher
        {
            unsigned operator()(const face& obj, unsigned seed=0) const
            {
                return obj.hash_code(seed);
            }
        };

        //- Symmetric hashing functor for face
        struct symmHasher
        {
            unsigned operator()(const face& obj, unsigned seed=0) const
            {
                return face::symmhash_code(obj, seed);
            }
        };


    // Housekeeping

        //- Identical to edge()
        Foam::edge faceEdge(label edgei) const { return this->edge(edgei); }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Hash face as contiguous (non-commutative) list data
template<> struct Hash<face> : face::hasher {};


//- Specialization to offset faces, used in ListListOps::combineOffset
template<>
struct offsetOp<face>
{
    face operator()(const face& x, const label offset) const
    {
        face f(x);
        f += offset;
        return f;
    }
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const face& a, const face& b);
inline bool operator!=(const face& a, const face& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceI.H"

#ifdef NoRepository
    #include "faceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
