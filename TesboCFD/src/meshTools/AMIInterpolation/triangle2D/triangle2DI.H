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

#include "OFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::triangle2D::writeOBJ
(
    Ostream& os,
    const triangle2D& tri,
    label offset
)
{
    os  << "v " << tri[0].x() << " " << tri[0].y() << " 0" << nl
        << "v " << tri[1].x() << " " << tri[1].y() << " 0" << nl
        << "v " << tri[2].x() << " " << tri[2].y() << " 0" << nl
        << "l"
        << " " << 1 + offset
        << " " << 2 + offset
        << " " << 3 + offset
        << " " << 1 + offset
        << endl;
}


inline Foam::label Foam::triangle2D::nClosePoints
(
    const triangle2D& triA,
    const triangle2D& triB
)
{
    label nClose = 0;

    FixedList<bool, 3> match(true);

    for (auto& a : triA)
    {
        forAll(triB, tb)
        {
            if (match[tb] && a.isClose(triB[tb]))
            {
                match[tb] = false;
                ++nClose;
                break;
            }
        }
    }

    return nClose;
}


inline Foam::scalar Foam::triangle2D::area
(
    const vector2D& a,
    const vector2D& b,
    const vector2D& c
)
{
    const vector2D e1(a - c);
    const vector2D e2(b - c);

    return 0.5*e1.perp(e2);
}


inline Foam::vector2D Foam::triangle2D::centre() const
{
    const triangle2D& tri = *this;

    return (1.0/3.0)*(tri[0] + tri[1] + tri[2]);
}


inline bool Foam::triangle2D::lineSegmentIntersectionPoint
(
    const vector2D& lp1,
    const vector2D& lp2,
    const vector2D& sp1,
    const vector2D& sp2,
    vector2D& intersection
)
{
    const vector2D r(lp2 - lp1);
    const vector2D s(sp2 - sp1);

    const scalar rcs = r.perp(s);

    if (mag(rcs) > ROOTVSMALL)
    {
        const scalar u = (sp1 - lp1).perp(r)/rcs;

        if (u >= -relTol && u <= 1+relTol)
        {
            intersection = sp1 + u*s;
            return true;
        }
    }

/*
        // Collapsed to line
        if (mag(triangle2D(lp1, lp2, sp1).area()) < absTol)
        {
            intersection = sp1;
            return true;
        }

        if (mag(triangle2D(lp1, lp2, sp2).area()) < absTol)
        {
            intersection = sp2;
            return true;
        }
*/

    if (debug)
    {
        OFstream os("bad-intersection.obj");
        os  << "v " << lp1.x() << " " << lp1.y() << " 0" << nl
            << "v " << lp2.x() << " " << lp2.y() << " 0" << nl
            << "v " << sp1.x() << " " << sp1.y() << " 0" << nl
            << "v " << sp2.x() << " " << sp2.y() << " 0" << nl
            << "l 1 2"
            << "l 3 4"
            << endl;
    }

    return false;
}


inline bool Foam::triangle2D::lineIntersectionPoint
(
    const vector2D& a,
    const vector2D& b,
    const vector2D& c,
    const vector2D& d,
    vector2D& intersection
)
{
    return lineSegmentIntersectionPoint(c, d, a, b, intersection);
}


inline bool Foam::triangle2D::lineIntersects
(
    const vector2D& a,
    const vector2D& b,
    const vector2D& c,
    const vector2D& d
)
{
    if
    (
        (triangle2D(a, c, d).order() == triangle2D(b, c, d).order())
     || (triangle2D(a, b, c).order() == triangle2D(a, b, d).order())
    )
    {
        return false;
    }


    DebugInfo<< "line " << a << b << " intersects " << c << d << endl;

    return true;
}


inline Foam::scalar Foam::triangle2D::area() const noexcept
{
    return area_;
}


inline Foam::label Foam::triangle2D::order() const
{
    return mag(area_) < SMALL ? 0 : sign(area_);
}


inline bool Foam::triangle2D::contains(const triangle2D& tri) const
{
    return
        pointInside(tri[0])
     && pointInside(tri[1])
     && pointInside(tri[2]);
}


inline bool Foam::triangle2D::isSame(const triangle2D& triB) const
{
    const triangle2D& triA = *this;

    return
        triA[0].isClose(triB[0])
     && triA[1].isClose(triB[1])
     && triA[2].isClose(triB[2]);
}


inline bool Foam::triangle2D::pointInside(const vector2D& p) const
{
    const triangle2D& triA = *this;

    for (int i = 0; i < 3; ++i)
    {
        if ((triA[(i + 1) % 3] - triA[i]).perp(p - triA[i]) < 0)
        {
            return false;
        }
    }

    return true;
}


// ************************************************************************* //
