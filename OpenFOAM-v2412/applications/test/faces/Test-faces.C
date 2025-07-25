/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

Application
    Test-faces

Description
    Simple tests for various faces

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelledTri.H"
#include "edge.H"
#include "faceList.H"
#include "triFaceList.H"
#include "pointList.H"
#include "ListOps.H"

using namespace Foam;


template<class Face>
void faceInfo(const Face& f, const UList<point>& points)
{
    Info<< f
        << " points:" << f.points(points)
        << " normal:" << f.unitNormal(points);
}


template<class Face>
void testSign
(
    const Face& f,
    const UList<point>& points,
    const UList<point>& testPoints
)
{
    for (const point& p : testPoints)
    {
        Info<< "  point:" << p << " sign=" << f.sign(p, points) << nl;
    }
}


template<class Face>
void testEdges(const Face& f)
{
    const label nEdges = f.nEdges();
    Info<< "face: " << f << nl
        << "flip: " << f.reverseFace() << nl
        << "  fc edges:" << flatOutput(f.edges()) << nl
        << "  rc edges:" << flatOutput(f.rcEdges()) << nl;

    Info<< "  forward edges" << nl;
    for (label edgei = 0; edgei < nEdges; ++edgei)
    {
        Info<< "    " << edgei << " : " << f.edge(edgei) << nl;
    }
    Info<< "  reverse edges" << nl;
    for (label edgei = 0; edgei < nEdges; ++edgei)
    {
        Info<< "    " << edgei << " : " << f.rcEdge(edgei) << nl;
    }
}


void testCompare(const triFace& a, const triFace& b)
{
    Info<< "compare: " << a << " with " << b
        <<  " == " << triFace::compare(a, b) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    pointList points1
    ({
        { 0, 0, 0},
        { -1, -1, 1},
        { 1, -1, -1},
        { 1, 1, -1},
        { -1, 1, 1}
    });

    pointList points2 = ListOps::create<point>
    (
        points1,
        [](const point& p){ return point(p.x(), p.y(), -p.z()); }
    );

    pointList testPoints
    ({
        { -2, -2, -2},
        { -2, -2, 2},
        { 0, 0, 0},
        { 2, 2, -2},
        { 2, 2, 2}
    });

    face f1({1, 2, 3, 4});
    Info<< "face:"; faceInfo(f1, points1); Info << nl;
    testSign(f1, points1, testPoints);

    Info<< "find within " << f1 << nl;
    for
    (
        const edge& e
      : edgeList
        ({
            {0, 1},
            {2, 1},
            {4, 1},
            f1.edge(1),
            {4, 3},
            {8, 1},
        }))
    {
        Info<< "    edge:" << e << " = " << f1.find(e)
            << " direction:" << f1.edgeDirection(e)
            << " contains:" << f1.contains(e) << nl;
    }

    Info<< "face:"; faceInfo(f1, points2); Info << nl;
    testSign(f1, points2, testPoints);
    Info<< nl;

    // Initializer list
    // triFace t1({1, 2, 3});

    // Component-wise
    {
        edge e1(3, 2, 1);  // Inadvertent sort!!!
        Info<< "edge:" << e1 << nl;
    }

    // Component-wise
    triFace t1(1, 2, 3);
    Info<< "triFace:"; faceInfo(t1, points1); Info << nl;
    testSign(t1, points1, testPoints);

    Info<< "find within " << t1 << nl;
    for
    (
        const edge& e
      : edgeList
        ({
            {0, 1},
            {2, 1},
            {4, 1},
            t1.edge(1),
            {4, 3},
            {8, 1},
        }))
    {
        Info<< "    edge:" << e << " = " << t1.find(e)
            << " direction:" << t1.edgeDirection(e)
            << " contains:" << t1.contains(e) << nl;
    }

    {
        scalarField fld({0, 20, 20, 30});
        Info<< "avg:" << t1.average(pointField::null(), fld) << nl;
    }

    Info<< "triFace:"; faceInfo(t1, points2); Info << nl;
    testSign(t1, points2, testPoints);
    Info<< nl;


    f1 = t1;
    Info<< "face:" << f1 << nl;

    f1 = t1.triFaceFace();
    Info<< "face:" << f1 << nl;

    #if 0
    // Expect failure, but triggers abort which cannot be caught
    const bool oldThrowingError = FatalError.throwing(true);
    try
    {
        labelledTri l1({1, 2, 3, 10, 24});
        Info<< "labelled:" << l1 << nl;
    }
    catch (const Foam::error& err)
    {
        WarningInFunction
            << "Caught FatalError " << err << nl << endl;
    }
    FatalError.throwing(oldThrowingError);
    #endif

    labelledTri l2({1, 2, 3});
    Info<< "labelled:" << l2 << nl;

    labelledTri l3({1, 2, 3, 10});
    Info<< "labelled:" << l3 << nl;

    t1.flip();
    l3.flip();

    Info<< "flip:" << t1 << nl;
    Info<< "flip:" << l3 << nl;

    {
        triFaceList faceList1
        ({
            triFace{1, 2, 3},
            triFace{4, 2, 100},
            triFace{1, 3, 2},
        });

        Info<< nl << "Test edges" << nl;

        for (const auto& f : faceList1)
        {
            testEdges(f);
            Info<< nl;
        }
    }

    {
        faceList faceList1
        ({
            face{1, 2, 3, 4},
            face{1, 4, 3, 2},
            face{4, 2, 100, 8, 35},
        });

        Info<< nl << "Test edges" << nl;

        for (const auto& f : faceList1)
        {
            testEdges(f);
            Info<< nl;
        }
    }

    {
        triFaceList faceList1
        ({
            triFace{1, 2, 3},
            triFace{1, 3, 2},
            triFace{3, 1, 2},
            triFace{4, 5, 1},
        });

        Info<< nl << "Test triFace compare" << nl;
        for (const triFace& a : faceList1)
        {
            for (const triFace& b : faceList1)
            {
                testCompare(a, b);
            }
        }
        Info<< nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
