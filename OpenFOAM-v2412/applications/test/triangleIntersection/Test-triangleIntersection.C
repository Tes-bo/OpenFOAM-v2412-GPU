/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Description
    Test bounding box / triangle intersection

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "line.H"
#include "Random.H"
#include "triangle.H"
#include "triangleFuncs.H"
#include "treeBoundBox.H"
#include "ListOps.H"

using namespace Foam;

//- simple helper to create a cube
boundBox cube(scalar start, scalar width)
{
    return boundBox
    (
        point::uniform(start),
        point::uniform(start + width)
    );
}


void printEdges(const treeBoundBox& bb)
{
    pointField pts(bb.points());

    for (const edge& e : treeBoundBox::edges)
    {
        Info<< pts[e.first()] << " -> " << pts[e.second()] << nl;
    }
}


void testIntersect(const treeBoundBox& bb, const triPoints& tri)
{
    int ninside = 0;

    if (bb.contains(tri.a())) ++ninside;
    if (bb.contains(tri.b())) ++ninside;
    if (bb.contains(tri.c())) ++ninside;

    Info<< "box: " << bb << endl;
    Info<< "tri: " << tri << endl;

    Info<< "num inside: " << ninside << nl;
    Info<< "intersects: " << bb.intersects(tri.tri())
        << ' ' << triangleFuncs::intersectBb(tri.tri(), bb) << nl << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    // Info<<"tree-bb faces: " << treeBoundBox::faces << nl
    //     <<"tree-bb edges: " << treeBoundBox::edges << endl;

    treeBoundBox bb(zero_one{});

    triPoints tri;
    tri.a() = point(-0.1, 0.5, 0.5);
    tri.b() = point(0.1, 0.6, 0.5);
    tri.c() = point(0.2, 0.4, 0.5);

    testIntersect(bb, tri);

    tri.a().z() = 1.1;
    tri.b().z() = 1.1;
    tri.c().z() = 1.1;

    testIntersect(bb, tri);

    tri.a().z() = 1 + 1e-15;
    tri.b().z() = 1 + 1e-15;
    tri.c().z() = 1 + 1e-15;

    testIntersect(bb, tri);

    tri.a() = point(-1, -1, -1);
    tri.b() = point(2, 2, -2);
    tri.c() = point(0, 0, 2);

    testIntersect(bb, tri);

    tri.a() = point(0.9, 1.1, 0);
    tri.b() = point(1.1, 0.9, 0);
    tri.c() = point(1, 1, 1.1);

    testIntersect(bb, tri);


    tri.a() = point(-1e-3, 1e-3,  -1);
    tri.b() = point( 1e-3, -1e-3, -1);
    tri.c() = point(0, 0, 2);

    testIntersect(bb, tri);

    tri.a() = point(-1e-3, 1e-3,  -1);
    tri.b() = point( 1e-3, -1e-3, -1);
    tri.c() = point(-1e-4, -1e-4, 2);

    testIntersect(bb, tri);

    return 0;
}


// ************************************************************************* //
