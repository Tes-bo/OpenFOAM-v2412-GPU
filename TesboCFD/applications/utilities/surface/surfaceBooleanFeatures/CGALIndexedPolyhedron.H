/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

Typedefs
    IndexedPolyhedron

Description
    CGAL data structures used for triSurface handling

    Define CGAL_INEXACT to use Exact_predicates_inexact_constructions kernel
    otherwise the more robust but much less efficient
    Exact_predicates_exact_constructions will be used.

\*---------------------------------------------------------------------------*/

#ifndef CGALIndexedPolyhedron_H
#define CGALIndexedPolyhedron_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include uint.H before CGAL headers to define __STDC_LIMIT_MACROS
#include "uint.H"

#include "CGAL3DKernel.H"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>

#include "label.H"


typedef CGAL::Point_3<K>  Point;
typedef CGAL::Segment_3<K> Segment;
typedef CGAL::Direction_3<K>  Direction;
typedef CGAL::Plane_3<K>  Plane;
typedef CGAL::Triangle_3<K>  Triangle;

// Define new class with color and define the polyhedron types
template<class Refs>
struct IndexedFace
:
    public CGAL::HalfedgeDS_face_base<Refs>
{
    Foam::label index;
    Foam::label region;
};

struct My_items
:
    public CGAL::Polyhedron_items_3
{
    template<class Refs, class Traits>
    struct Face_wrapper
    {
        typedef IndexedFace<Refs> Face;
    };
};


//typedef CGAL::Polyhedron_3<K>  Polyhedron;
typedef CGAL::Polyhedron_3<K, My_items>  Polyhedron;

typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef Polyhedron::Edge_iterator Edge_iterator;
typedef Polyhedron::Vertex Vertex;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Edge_iterator Edge_iterator;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator HFCC;
typedef Polyhedron::Vertex_const_iterator                  VCI;

typedef CGAL::Nef_polyhedron_3<K>  Nef_polyhedron;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
