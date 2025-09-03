/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::faceAreaIntersect::triangleFan
(
    const face& f,
    DynamicList<face>& faces
)
{
    if (f.size() > 2)
    {
        const label v0 = 0;

        face indices(3);

        for (label i = 1; i < f.size() - 1; ++i)
        {
            indices[0] = f[v0];
            indices[1] = f[i];
            indices[2] = f[i + 1];
            faces.append(indices);
        }
    }
}


inline void Foam::faceAreaIntersect::setTriPoints
(
    const point& a,
    const point& b,
    const point& c,
    label& count,
    FixedList<triPoints, 10>& tris
) const
{
    triPoints& tp = tris[count++];
    tp[0] = a;
    tp[1] = b;
    tp[2] = c;
}


inline Foam::triPoints Foam::faceAreaIntersect::getTriPoints
(
    const pointField& points,
    const face& f,
    const bool reverse
) const
{
    triPoints result;

    if (reverse)
    {
        result[2] = points[f[0]];
        result[1] = points[f[1]];
        result[0] = points[f[2]];
    }
    else
    {
        result[0] = points[f[0]];
        result[1] = points[f[1]];
        result[2] = points[f[2]];
    }

    return result;
}


inline Foam::point Foam::faceAreaIntersect::planeIntersection
(
    const FixedList<scalar, 3>& d,
    const triPoints& t,
    const label negI,
    const label posI
) const
{
    scalar dp = d[posI];
    scalar dn = d[negI];
    return (dp*t[negI] - dn*t[posI])/(-dn + dp);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar& Foam::faceAreaIntersect::tolerance()
{
    return tol;
}


bool Foam::faceAreaIntersect::cacheTriangulation() const
{
    return cacheTriangulation_;
}


const Foam::DynamicList<Foam::triPoints>
Foam::faceAreaIntersect::triangles() const
{
    return triangles_;
}


Foam::DynamicList<Foam::triPoints>& Foam::faceAreaIntersect::triangles()
{
    return triangles_;
}


// ************************************************************************* //
