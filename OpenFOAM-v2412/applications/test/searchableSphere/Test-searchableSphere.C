/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application
    Test-searchableSphere

Description
    Basic tests for searchable sphere

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "searchableSphere.H"
#include "unitConversion.H"
#include "Random.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

//- Evaluate using implicit form of the spheroid equation.
scalar evaluate(const point& p, const searchableSphere& sph)
{
    return
    (
        sqr((p.x() - sph.centre().x()) / sph.radii().x())
      + sqr((p.y() - sph.centre().y()) / sph.radii().y())
      + sqr((p.z() - sph.centre().z()) / sph.radii().z())
    );
}


void doTest1(const searchableSphere& sph)
{
    Info<< nl << "origin:" << sph.centre() << " radius:" << sph.radius();

    if (sph.shape() == searchableSphere::SPHERE)
    {
        Info<< " [sphere]" << nl;
    }
    else
    {
        Info<< " radii:" << sph.radii() << nl;
    }

//     boundBox bb(point(0, 0, 0), point(30, 30, 30));
//     Info<< "overlaps: " << Switch(sph.overlaps(bb)) << endl;

    Pair<scalar> angles(-pi/2, pi/2);

    point surfPt = sph.surfacePoint(angles.first(), angles.second());
    vector norm = sph.surfaceNormal(angles.first(), angles.second());

    Info<< "point at ("
        << radToDeg(angles.first()) << ' '
        << radToDeg(angles.second()) << ") deg" << nl;

    Info<< "surface" << nl
        << "    eval: " << evaluate(surfPt, sph) << nl
        << "   point:" << surfPt << nl
        << "  normal:" << norm << nl;

    {
        List<pointIndexHit> hits(1);
        vectorField norms;

        hits[0].hitPoint(surfPt, 0);

        sph.getNormal(hits, norms);

        Info<< "  normal:" << norms[0] << nl;
    }


    Random rndGen;

    point testPt1 =
    (
        surfPt + 0.1*sph.radius()*norm
      + 0.01*rndGen.sample01<vector>()
    );


    // Scale by max radius and shift by origin
    const auto adjustPoint =
        [&](point& p) -> void
        {
            p *= sph.radius();
            p += sph.centre();
        };


    List<pointIndexHit> hits;
    pointField query
    ({
        testPt1,
        point(-2, -2, -2)
    });

    forAll(query, pointi)
    {
        if (pointi) adjustPoint(query[pointi]);
    }

    sph.findNearest
    (
        query,
        scalarField(query.size(), GREAT),
        hits
    );

    Info<< "query:" << nl;

    forAll(hits, pointi)
    {
        // Should never miss
        Info<< "  "
            << query[pointi] << " => "
            << hits[pointi].hitPoint() << nl;
    }


    pointField lineBeg
    ({
        point(-2, -2, -2),
        point(0,0,0),
        point(0,0,0),
        point(0,0,0)
    });
    pointField lineEnd
    ({
        point(2, 2, 2),
        point(2, 0, 0),
        point(0, 2, 0),
        point(0, 0, 4)
    });

    for (point& p : lineBeg)
    {
        adjustPoint(p);
    }

    for (point& p : lineEnd)
    {
        adjustPoint(p);
    }

    sph.findLineAny
    (
        lineBeg,
        lineEnd,
        hits
    );

    Info<< "lines:" << nl;

    forAll(hits, pointi)
    {
        // Should never miss
        Info<< "  "
            << lineBeg[pointi] << " -> "
            << lineEnd[pointi] << " = "
            << hits[pointi].hitPoint() << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::noBanner();
    #include "setRootCase.H"

    // Use dummy Time for objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New());

    const IOobject io
    (
        "sphere",
        *dummyTimePtr,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    Info<< "Testing searchable sphere" << endl;

    doTest1
    (
        searchableSphere
        (
            io,
            point(0.5, 0.5, 0.5),
            scalar(2)
        )
    );

    doTest1
    (
        searchableSphere
        (
            io,
            point(0.5, 0.5, 0.5),
            vector(1.999, 2, 2.001)
        )
    );

    doTest1
    (
        searchableSphere
        (
            io,
            point(0, 0, 0),
            vector(3, 3, 4)
        )
    );

    Info<< "\nDone\n" << endl;
    return 0;
}

// ************************************************************************* //
