/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::edge

Description
    An edge is a list of two vertex labels.
    This can correspond to a directed graph edge or an edge on a mesh.

    The edge is implemented as a Pair/FixedList of labels.
    As well as geometrically relevant methods, it also provides methods
    similar to HashSet for additional convenience.
    Valid vertex labels are always non-negative (eg, since they correspond to
    addressing within the mesh). The value '-1' is used to tag invalid
    point labels that correspond conceptually to open 'slots', which
    can be filled with a HashSet-like functionality.

SourceFiles
    edge.C
    edgeI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_edge_H
#define Foam_edge_H

#include "labelPair.H"
#include "line.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class edge Declaration
\*---------------------------------------------------------------------------*/

class edge
:
    public labelPair
{
public:

    // Static Data Members

        //- The typeName ("edge")
        static const char* const typeName;


    // Constructors

        //- Default construct, with invalid vertex labels (-1)
        inline edge();

        //- Construct from two vertex labels
        inline edge(const label from, const label to);

        //- Construct from pair of vertex labels
        inline edge(const labelPair& pair);

        //- Construct from list of vertex labels
        inline edge(const FixedList<label, 2>& list);

        //- Construct from two vertex labels, sorted with first less-than second
        inline edge(const label from, const label to, const bool doSort);

        //- Construct from list, sorted with first less-than second
        inline edge(const FixedList<label, 2>& list, const bool doSort);

        //- Copy construct from a subset of vertex labels
        inline edge
        (
            const UList<label>& list,
            const FixedList<label, 2>& indices
        );

        //- Construct from Istream
        inline edge(Istream& is);


    // Member Functions

    // Access

        //- The first vertex
        label a() const noexcept { return labelPair::first(); }

        //- The second vertex
        label b() const noexcept { return labelPair::second(); }

        //- The first vertex
        label& a() noexcept { return labelPair::first(); }

        //- The second vertex
        label& b() noexcept { return labelPair::second(); }

        //- The start (first) vertex label
        label start() const noexcept { return labelPair::first(); }

        //- The end (last/second) vertex label
        label end() const noexcept { return labelPair::second(); }

        //- The start (first) vertex label
        label& start() noexcept { return labelPair::first(); }

        //- The end (second/last) vertex label
        label& end() noexcept { return labelPair::second(); }


        //- Return reverse edge as copy.
        //  No special handling of negative vertex labels.
        inline edge reverseEdge() const;


    // Queries

        //- Return the smallest vertex label used by the edge
        //  No special handling of negative vertex labels.
        inline label minVertex() const noexcept;

        //- Return the largest vertex label used by the edge
        //  No special handling of negative vertex labels.
        inline label maxVertex() const noexcept;

        //- True if the vertices are unique and non-negative.
        inline bool good() const noexcept;

        //- Return true if the vertex label is contained in the edge.
        //  Always false for a negative vertex label.
        inline bool contains(const label vertex) const noexcept;

        //- Return local index (0,1) of vertex label in edge -1 on failure
        //  Always return -1 for a negative vertex label.
        inline label which(const label vertex) const;

        //- True if the edge has at least one vertex in common with other
        inline bool connected(const edge& other) const;

        //- Return vertex common with other edge or -1 on failure
        //  Negative vertex labels are never considered common between edges.
        inline label commonVertex(const edge& other) const;

        //- Given one vertex label, return the other one.
        //  No special treatment for negative vertex labels.
        inline label otherVertex(const label vertex) const;

        //- Do the edges share a common vertex index?
        //  Negative vertex labels never connect.
        bool connects(const edge& other) const { return connected(other); }


    // Editing

        //- 'Collapse' edge by marking duplicate vertex labels as '-1',
        //- the lower vertex is retained.
        //  Return the effective size after collapsing.
        inline label collapse();


    // Hash-like Functions

        //- Return the number of unique, valid (non -1) vertex labels.
        //  Similar to a HashTable::size().
        inline label count() const;

        //- Return true if edge has no valid vertex labels.
        inline bool empty() const noexcept;

        //- 'Clears' edge by setting both ends to invalid vertex labels.
        inline void clear();


        //- Fill any open slot with the vertex label
        //- (if not previously contained in the edge).
        //  Returns true on success. A negative vertex label never inserts.
        //  Similar to a HashTable::insert().
        inline bool insert(const label vertex);

        //- Insert values, using begin/end iterators.
        template<class InputIterator>
        inline label insert(InputIterator begIter, InputIterator endIter);

        //- Fill open slots with the indices if they did not previously exist.
        //  Returns true on success. Negative labels never insert.
        //  Return the number of slots filled.
        //  Similar to a HashTable::insert().
        inline label insert(std::initializer_list<label> list);

        //- Fill open slots with the indices if they did not previously exist.
        //  Returns true on success. Negative labels never insert.
        //  Return the number of slots filled.
        //  Similar to a HashTable::insert().
        template<unsigned N>
        inline label insert(const FixedList<label, N>& list);

        //- Fill open slots with the indices if they did not previously exist.
        //  Returns true on success. Negative labels never insert.
        //  Return the number of slots filled.
        //  Similar to a HashTable::insert().
        inline label insert(const labelUList& list);


        //- Remove an existing vertex from the edge and set its location
        //- to '-1'. A negative vertex label never removes.
        //  Returns the number of changes.
        //  Similar to a HashTable::erase().
        inline label erase(const label vertex);

        //- Remove values, using begin/end iterators.
        template<class InputIterator>
        inline label erase(InputIterator begIter, InputIterator endIter);

        //- Remove existing indices from the edge and set locations to '-1'.
        //  Returns the number of changes.
        inline label erase(std::initializer_list<label> list);

        //- Remove existing indices from the edge and set locations to '-1'.
        //  Returns the number of changes.
        template<unsigned N>
        inline label erase(const FixedList<label, N>& list);

        //- Remove existing indices from the edge and set locations to '-1'.
        //  Returns the number of changes.
        inline label erase(const labelUList& list);


    // Geometric Functions

        //- Return centre point (centroid) of the edge.
        //  No special handling of negative vertex labels.
        inline point centre(const UList<point>& pts) const;

        //- Return the vector (from first to second).
        //  No special handling of negative vertex labels.
        inline vector vec(const UList<point>& pts) const;

        //- Return the unit vector (from first to second).
        //  No special handling of negative vertex labels.
        inline vector unitVec(const UList<point>& pts) const;

        //- The length (L2-norm) of the edge vector.
        //  No special handling of negative vertex labels.
        inline scalar mag(const UList<point>& pts) const;

        //- The length (L2-norm) squared of the edge vector.
        //  No special handling of negative vertex labels.
        inline scalar magSqr(const UList<point>& pts) const;

        //- The enclosing (bounding) box for the edge
        inline Pair<point> box(const UList<point>& pts) const;

        //- Return edge line
        //  No special handling of negative vertex labels.
        inline linePointRef line(const UList<point>& pts) const;


    // Comparison

        //- Compare edges
        //  \return
        //  -  0: different
        //  - +1: identical values and order used
        //  - -1: identical values, but in different order
        static inline int compare(const edge& a, const edge& b);


    // Member Operators

        //- Return edge element. Index should be limited to 0/1.
        inline label& operator[](const label i);

        //- Return constant edge element. Index should be limited to 0/1.
        inline const label& operator[](const label i) const;


    // Hashing

        //- The (commutative) hash value for edge, hashes lower value first
        inline unsigned hash_code(unsigned seed=0) const
        {
            Foam::Hash<label> op;
            if (second() < first())
            {
                return op(first(), op(second(), seed));
            }
            else
            {
                return op(second(), op(first(), seed));
            }
        }

        //- Hashing functor for edge (commutative)
        //  Also useful for inheritance in sub-classes
        struct hasher
        {
            unsigned operator()(const edge& obj, unsigned seed=0) const
            {
                return obj.hash_code(seed);
            }
        };


    // Housekeeping

        //- Same as contains()
        bool found(label vertex) const { return contains(vertex); }

        //- Same as good()
        bool valid() const noexcept { return good(); }

        //- Deprecated(2021-04) hashing functor. Use hasher()
        // \deprecated(2021-04) - use hasher() functor
        template<class Unused=bool>
        struct Hash : edge::hasher
        {
            FOAM_DEPRECATED_FOR(2021-04, "hasher()") Hash() {}
        };
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for edge (a pair of labels)
template<> struct is_contiguous<edge> : std::true_type {};

//- Contiguous label data for edge (a pair of labels)
template<> struct is_contiguous_label<edge> : std::true_type {};

//- Hashing for edge uses commutative (incremental) hash
template<> struct Hash<edge> : edge::hasher {};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return reverse of an edge
inline edge reverse(const edge& e)
{
    return edge(e.second(), e.first());
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

//- Compare edges for equal content, ignoring orientation
inline bool operator==(const edge& a, const edge& b);

//- Compare edges for non-equal content, ignoring orientation
inline bool operator!=(const edge& a, const edge& b);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
