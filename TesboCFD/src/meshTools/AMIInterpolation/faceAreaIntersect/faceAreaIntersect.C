/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "faceAreaIntersect.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::faceAreaIntersect::triangulationMode
>
Foam::faceAreaIntersect::triangulationModeNames_
({
    { triangulationMode::tmFan, "fan" },
    { triangulationMode::tmMesh, "mesh" },
});

Foam::scalar Foam::faceAreaIntersect::tol = 1e-6;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::faceAreaIntersect::triSliceWithPlane
(
    const triPoints& tri,
    const plane& pln,
    FixedList<triPoints, 10>& tris,
    label& nTris,
    const scalar len
) const
{
    // distance to cutting plane
    FixedList<scalar, 3> d;

    // determine how many of the points are above the cutting plane
    label nCoPlanar = 0;
    label nPos = 0;
    label posI = -1;
    label negI = -1;
    label copI = -1;
    forAll(tri, i)
    {
        d[i] = pln.signedDistance(tri[i]);

        if (mag(d[i]) < tol*len)
        {
            nCoPlanar++;
            copI = i;
            d[i] = 0.0;
        }
        else
        {
            if (d[i] > 0)
            {
                nPos++;
                posI = i;
            }
            else
            {
                negI = i;
            }
        }
    }


    // Determine triangle area contribution

    if
    (
        (nPos == 3)
     || ((nPos == 2) && (nCoPlanar == 1))
     || ((nPos == 1) && (nCoPlanar == 2))
    )
    {
        /*
                /\          _____
               /  \         \   /          /\
              /____\         \ /          /  \
            __________    ____v____    __/____\__

            all points above cutting plane
            - add complete triangle to list
        */

        tris[nTris++] = tri;
    }
    else if ((nPos == 2) && (nCoPlanar == 0))
    {
        /*
            i1________i2
              \      /
             --\----/--
                \  /
                 \/
                 i0

            2 points above plane, 1 below
            - resulting quad above plane split into 2 triangles
            - forget triangle below plane
        */

        // point under the plane
        label i0 = negI;

        // indices of remaining points
        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);

        // determine the two intersection points
        point p01 = planeIntersection(d, tri, i0, i1);
        point p02 = planeIntersection(d, tri, i0, i2);

        // forget triangle below plane
        // - decompose quad above plane into 2 triangles and add to list
        setTriPoints(tri[i1], tri[i2], p02, nTris, tris);
        setTriPoints(tri[i1], p02, p01, nTris, tris);
    }
    else if (nPos == 1)
    {
        // point above the plane
        label i0 = posI;

        if (nCoPlanar == 0)
        {
            /*
                     i0
                     /\
                    /  \
                 --/----\--
                  /______\
                i2        i1

                1 point above plane, 2 below
                - keep triangle above intersection plane
                - forget quad below plane
            */

            // indices of remaining points
            label i1 = d.fcIndex(i0);
            label i2 = d.fcIndex(i1);

            // determine the two intersection points
            point p01 = planeIntersection(d, tri, i1, i0);
            point p02 = planeIntersection(d, tri, i2, i0);

            // add triangle above plane to list
            setTriPoints(tri[i0], p01, p02, nTris, tris);
        }
        else
        {
            /*
                  i0
                  |\
                  | \
                __|__\_i2_
                  |  /
                  | /
                  |/
                  i1

                1 point above plane, 1 on plane, 1 below
                - keep triangle above intersection plane
            */

            // point indices
            label i1 = negI;
            label i2 = copI;

            // determine the intersection point
            point p01 = planeIntersection(d, tri, i1, i0);

            // add triangle above plane to list - clockwise points
            if (d.fcIndex(i0) == i1)
            {
                setTriPoints(tri[i0], p01, tri[i2], nTris, tris);
            }
            else
            {
                setTriPoints(tri[i0], tri[i2], p01, nTris, tris);
            }
        }
    }
    else
    {
        /*
            _________    __________    ___________
                             /\          \    /
               /\           /  \          \  /
              /  \         /____\          \/
             /____\

            all points below cutting plane - forget
        */
    }
}


void Foam::faceAreaIntersect::triangleIntersect
(
    const triPoints& src,
    const point& tgt0,
    const point& tgt1,
    const point& tgt2,
    const vector& n,
    scalar& area,
    vector& centroid
) const
{
    // Work storage
    FixedList<triPoints, 10> workTris1;
    label nWorkTris1 = 0;

    FixedList<triPoints, 10> workTris2;
    label nWorkTris2 = 0;

    // Cut source triangle with all inward pointing faces of target triangle
    // - triangles in workTris1 are inside target triangle

    const scalar srcArea(src.mag());
    if (srcArea < ROOTVSMALL)
    {
        return;
    }

    // Typical length scale
    const scalar t = sqrt(srcArea);

    // Edge 0
    {
        // Cut triangle src with plane and put resulting sub-triangles in
        // workTris1 list

        scalar s = mag(tgt1 - tgt0);
        if (s < ROOTVSMALL)
        {
            return;
        }

        // Note: outer product with n pre-scaled with edge length. This is
        //       purely to avoid numerical errors because of drastically
        //       different vector lengths.
        const vector n0((tgt0 - tgt1)^(-s*n));
        const scalar magSqrN0(magSqr(n0));
        if (magSqrN0 < ROOTVSMALL)
        {
            // Triangle either zero edge length (s == 0) or
            // perpendicular to face normal n. In either case zero
            // overlap area
            return;
        }
        else
        {
            plane pl0(tgt0, n0/Foam::sqrt(magSqrN0), false);
            triSliceWithPlane(src, pl0, workTris1, nWorkTris1, t);
        }
    }

    if (nWorkTris1 == 0)
    {
        return;
    }

    // Edge 1
    {
        // Cut workTris1 with plane and put resulting sub-triangles in
        // workTris2 list (re-use tris storage)

        scalar s = mag(tgt2 - tgt1);
        if (s < ROOTVSMALL)
        {
            return;
        }
        const vector n1((tgt1 - tgt2)^(-s*n));
        const scalar magSqrN1(magSqr(n1));

        if (magSqrN1 < ROOTVSMALL)
        {
            // Triangle either zero edge length (s == 0) or
            // perpendicular to face normal n. In either case zero
            // overlap area
            return;
        }
        else
        {
            plane pl1(tgt1, n1/Foam::sqrt(magSqrN1), false);

            nWorkTris2 = 0;

            for (label i = 0; i < nWorkTris1; ++i)
            {
                triSliceWithPlane(workTris1[i], pl1, workTris2, nWorkTris2, t);
            }

            if (nWorkTris2 == 0)
            {
                return;
            }
        }
    }

    // Edge 2
    {
        // Cut workTris2 with plane and put resulting sub-triangles in
        // workTris1 list (re-use workTris1 storage)

        scalar s = mag(tgt2 - tgt0);
        if (s < ROOTVSMALL)
        {
            return;
        }
        const vector n2((tgt2 - tgt0)^(-s*n));
        const scalar magSqrN2(magSqr(n2));

        if (magSqrN2 < ROOTVSMALL)
        {
            // Triangle either zero edge length (s == 0) or
            // perpendicular to face normal n. In either case zero
            // overlap area
            return;
        }
        else
        {
            plane pl2(tgt2, n2/Foam::sqrt(magSqrN2), false);

            nWorkTris1 = 0;

            for (label i = 0; i < nWorkTris2; ++i)
            {
                triSliceWithPlane(workTris2[i], pl2, workTris1, nWorkTris1, t);
            }

            if (nWorkTris1 == 0)
            {
                return;
            }
            else
            {
                // Calculate area of sub-triangles
                for (label i = 0; i < nWorkTris1; ++i)
                {
                    // Area of intersection
                    const scalar currArea = workTris1[i].mag();
                    area += currArea;

                    // Area-weighted centroid of intersection
                    centroid += currArea*workTris1[i].centre();

                    if (cacheTriangulation_)
                    {
                        triangles_.append(workTris1[i]);
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::faceAreaIntersect::faceAreaIntersect
(
    const pointField& pointsA,
    const pointField& pointsB,
    const DynamicList<face>& trisA,
    const DynamicList<face>& trisB,
    const bool reverseB,
    const bool cacheTriangulation
)
:
    pointsA_(pointsA),
    pointsB_(pointsB),
    trisA_(trisA),
    trisB_(trisB),
    reverseB_(reverseB),
    cacheTriangulation_(cacheTriangulation),
    triangles_(cacheTriangulation ? 10 : 0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceAreaIntersect::triangulate
(
    const face& f,
    const pointField& points,
    const triangulationMode& triMode,
    faceList& faceTris
)
{
    faceTris.resize(f.nTriangles());

    switch (triMode)
    {
        case triangulationMode::tmFan:
        {
            for (label i = 0; i < f.nTriangles(); ++i)
            {
                faceTris[i] = face(3);
                faceTris[i][0] = f[0];
                faceTris[i][1] = f[i + 1];
                faceTris[i][2] = f[i + 2];
            }

            break;
        }
        case triangulationMode::tmMesh:
        {
            const label nFaceTris = f.nTriangles();

            label nFaceTris1 = 0;
            const label nFaceTris2 = f.triangles(points, nFaceTris1, faceTris);

            if (nFaceTris != nFaceTris1 || nFaceTris != nFaceTris2)
            {
                FatalErrorInFunction
                    << "The numbers of reported triangles in the face do not "
                    << "match that generated by the triangulation"
                    << exit(FatalError);
            }
        }
    }
}


void Foam::faceAreaIntersect::calc
(
    const face& faceA,
    const face& faceB,
    const vector& n,
    scalar& area,
    vector& centroid
) const
{
    if (cacheTriangulation_)
    {
        triangles_.clear();
    }

    area = 0.0;
    centroid = vector::zero;

    // Intersect triangles
    for (const face& triA : trisA_)
    {
        triPoints tpA = getTriPoints(pointsA_, triA, false);

        for (const face& triB : trisB_)
        {
            if (reverseB_)
            {
                triangleIntersect
                (
                    tpA,
                    pointsB_[triB[0]],
                    pointsB_[triB[1]],
                    pointsB_[triB[2]],
                    n,
                    area,
                    centroid
                );
            }
            else
            {
                triangleIntersect
                (
                    tpA,
                    pointsB_[triB[2]],
                    pointsB_[triB[1]],
                    pointsB_[triB[0]],
                    n,
                    area,
                    centroid
                );
            }
        }
    }

    // Area weighed centroid
    if (area > 0)
    {
        centroid /= area;
    }
}


bool Foam::faceAreaIntersect::overlaps
(
    const face& faceA,
    const face& faceB,
    const vector& n,
    const scalar threshold
) const
{
    scalar area = 0.0;
    vector centroid(Zero);

    // Intersect triangles
    for (const face& triA : trisA_)
    {
        const triPoints tpA = getTriPoints(pointsA_, triA, false);

        for (const face& triB : trisB_)
        {
            if (reverseB_)
            {
                triangleIntersect
                (
                    tpA,
                    pointsB_[triB[0]],
                    pointsB_[triB[1]],
                    pointsB_[triB[2]],
                    n,
                    area,
                    centroid
                );
            }
            else
            {
                triangleIntersect
                (
                    tpA,
                    pointsB_[triB[2]],
                    pointsB_[triB[1]],
                    pointsB_[triB[0]],
                    n,
                    area,
                    centroid
                );
            }

            if (area > threshold)
            {
                return true;
            }
        }
    }

    return false;
}


// ************************************************************************* //
