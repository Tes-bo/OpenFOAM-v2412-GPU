/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "PolyhedronReader.H"

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class HDS>
void Foam::PolyhedronReader::Build_triangle<HDS>::operator()(HDS& hds)
{
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

    B.begin_surface(s_.nPoints(), s_.size());

    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;

    for (const auto& pt : s_.points())
    {
        B.add_vertex(Point(pt.x(), pt.y(), pt.z()));
    }

    for (const auto& f : s_)
    {
        B.begin_facet();

        for (const label verti : f)
        {
            B.add_vertex_to_facet(verti);
        }

        B.end_facet();
    }

    B.end_surface();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PolyhedronReader::PolyhedronReader(const triSurface& s, Polyhedron& p)
{
    Build_triangle<HalfedgeDS> triangle(s);
    p.delegate(triangle);

    // Populate index and region
    Foam::label nTris = 0;

    for
    (
        Facet_iterator fi = p.facets_begin();
        fi != p.facets_end();
        ++fi
    )
    {
        fi->index = nTris;
        fi->region = s[nTris].region();

        ++nTris;
    }
}


// ************************************************************************* //
