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

\*---------------------------------------------------------------------------*/

#include "triangle2D.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::triangle2D::debug = 0;

Foam::scalar Foam::triangle2D::relTol = 1e-8;

Foam::scalar Foam::triangle2D::absTol = 1e-10;

Foam::FixedList<Foam::vector2D, 7> Foam::triangle2D::work_
(
    vector2D::zero
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triangle2D::triangle2D
(
    const vector2D& a,
    const vector2D& b,
    const vector2D& c,
    const bool orient
)
:
    FixedList<vector2D, 3>({a, b, c}),
    area_(area(a, b, c))
{
    if (orient && area_ < 0)
    {
        // Inverted triangle
        triangle2D& tri = *this;
        vector2D tmp = tri[0];
        tri[0] = tri[2];
        tri[2] = tmp;

        area_ = mag(area_);
    }
}


Foam::triangle2D::triangle2D
(
    const vector& a3d,
    const vector& b3d,
    const vector& c3d,
    const vector& axis1,
    const vector& axis2,
    const bool orient
)
:
    triangle2D
    (
        vector2D(axis1 & a3d, axis2 & a3d),
        vector2D(axis1 & b3d, axis2 & b3d),
        vector2D(axis1 & c3d, axis2 & c3d),
        orient
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::triangle2D::snapClosePoints(const triangle2D& triB)
{
    label nSnapped = 0;
    triangle2D& triA = *this;

    FixedList<bool, 3> match(true);

    for (auto& a : triA)
    {
        forAll(triB, tb)
        {
            if (match[tb] && a.isClose(triB[tb]))
            {
                a = triB[tb];
                match[tb] = false;
                ++nSnapped;
                break;
            }
        }
    }

    return nSnapped;
}


void Foam::triangle2D::interArea
(
    const triangle2D& triB,
    vector2D& centre,
    scalar& area
) const
{
    const triangle2D& triA = *this;

    // Potential short-cut if the triangles are the same-ish.  Happens rarely
    // for moving mesh cases.
    // if (nClosePoints(triA, triB) == 3)
    // {
    //     centre = triA.centre();
    //     area = triA.area();
    //     return;
    // }

    if (debug)
    {
        static int nInter = 0;
        OFstream os("intersection-tris-"+Foam::name(nInter++)+".obj");
        writeOBJ(os, triA, 0);
        writeOBJ(os, triB, 3);
        Info<< "written " << os.name() << endl;
    }


    // Use work_ to store intersections

    // Current number of intersections
    int nPoint = 0;

    static FixedList<vector2D, 7> work2;

    // Clipped triangle starts as triA
    work2[0] = triA[0];
    work2[1] = triA[1];
    work2[2] = triA[2];

    int nPoint2 = 3;


    vector2D intersection(Zero);

    // Cut triA using triB's edges as clipping planes
    for (int i0 = 0; i0 <= 2; ++i0)
    {
        if (debug)
        {
            Info<< "\nstarting points:";
            for (int i = 0; i < nPoint2; ++i)
            {
                Info<< work2[i];
            }
            Info<< endl;
        }

        // Clipping plane points
        const label i1 = (i0 + 1) % 3;
        const vector2D& c0 = triB[i0];
        const vector2D& c1 = triB[i1];

        nPoint = 0;

        // Test against all intersection poly points
        for (int j = 0; j < nPoint2; ++j)
        {
            const vector2D& p0 = work2[j];
            const vector2D& p1 = work2[(j+1) % nPoint2];

            if (triangle2D(c0, c1, p0).order() == 1)
            {
                if (triangle2D(c0, c1, p1).order() == 1)
                {
                    work_[nPoint++] = p1;
                }
                else
                {
                    if (lineIntersectionPoint(p0, p1, c0, c1, intersection))
                    {
                        work_[nPoint++] = intersection;
                    }
                }
            }
            else
            {
                if (triangle2D(c0, c1, p1).order() == 1)
                {
                    if (lineIntersectionPoint(p0, p1, c0, c1, intersection))
                    {
                        work_[nPoint++] = intersection;
                    }
                    work_[nPoint++] = p1;
                }
            }
        }

        work2 = work_;
        nPoint2 = nPoint;
    }

    if (debug)
    {
        static int n = 0;
        OFstream os("intersection-poly-"+Foam::name(n++)+".obj");
        for (int i = 0; i < nPoint; ++i)
        {
            os << "v " << work_[i].x() << " " << work_[i].y() << " 0" << nl;
        }
        os  << "l";
        for (int i = 0; i < nPoint; ++i)
        {
            os << " " << (i + 1);
        }
        os  << " 1" << endl;
        Info<< "written " << os.name() << endl;

        Info<< "Intersection points:" << endl;
        for (int i = 0; i < nPoint; ++i)
        {
            Info<< "    " << work_[i] << endl;
        }
    }

    // Calculate the area of the clipped triangle
    scalar twoArea = 0;
    centre = vector2D::zero;
    if (nPoint)
    {
        for (int i = 0; i < nPoint - 1; ++i)
        {
            twoArea += work_[i].x()*work_[i+1].y();
            twoArea -= work_[i].y()*work_[i+1].x();
            centre += work_[i];
        }
        twoArea += work_[nPoint-1].x()*work_[0].y();
        twoArea -= work_[nPoint-1].y()*work_[0].x();

        centre += work_[nPoint - 1];
        centre /= scalar(nPoint);
    }

    area = 0.5*twoArea;
}


Foam::scalar Foam::triangle2D::interArea(const triangle2D& triB) const
{
    vector2D dummyCentre(Zero);
    scalar area;

    interArea(triB, dummyCentre, area);

    return area;
}


bool Foam::triangle2D::overlaps(const triangle2D& triB) const
{
    const triangle2D& triA = *this;

    // True if any of triB's edges intersect a triA edge
    for (int i = 0; i < 3; ++i)
    {
        int i1 = (i + 1) % 3;

        for (int j = 0; j < 3; ++j)
        {
            int j1 = (j + 1) % 3;

            if (lineIntersects(triA[i], triA[i1], triB[j], triB[j1]))
            {
                return true;
            }
        }
    }

    return
        (nClosePoints(triA, triB) == 3) // same tri
     || triA.contains(triB)             // triA contains triB
     || triB.contains(triA);            // triB contains triA
}


// ************************************************************************* //
