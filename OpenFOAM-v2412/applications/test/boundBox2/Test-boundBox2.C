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
    Test bounding box behaviour

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "line.H"
#include "Random.H"
#include "treeBoundBox.H"
#include "bitSet.H"
#include "HashSet.H"
#include "ListOps.H"

using namespace Foam;

//- simple helper to create a cube, given lower corner and width
boundBox cube(scalar start, scalar width)
{
    return boundBox
    (
        point::uniform(start),
        point::uniform(start + width)
    );
}

//- simple helper to create a cube, given mid-point and width
boundBox cubeAt(const point& mid, scalar width)
{
    boundBox bb(mid);
    bb.grow(0.5*width);

    return bb;
}


word faceName(direction whichFace)
{
    switch (whichFace)
    {
        case treeBoundBox::LEFT   : return "-x";
        case treeBoundBox::RIGHT  : return "+x";

        case treeBoundBox::BOTTOM : return "-y";
        case treeBoundBox::TOP    : return "+y";

        case treeBoundBox::BACK   : return "-z";
        case treeBoundBox::FRONT  : return "+z";
    }

    return "??";
}


word octantName(direction octant)
{
    word str("-x-y-z");

    if (octant & treeBoundBox::RIGHTHALF)
    {
        str[0] = '+';
    }
    if (octant & treeBoundBox::TOPHALF)
    {
        str[2] = '+';
    }
    if (octant & treeBoundBox::FRONTHALF)
    {
        str[4] = '+';
    }
    return str;
}


void testOverlaps(const treeBoundBox& bb, const treeBoundBox& searchBox)
{
    FixedList<bool, 8> overlaps;

    for (direction octant = 0; octant < 8; ++octant)
    {
        overlaps[octant] = bb.subOverlaps(octant, searchBox);
    }

    Info<< "box " << bb << " and " << searchBox << nl;

    Info<< "overlaps any:" << bb.overlaps(searchBox)
        << " octants: " << overlaps << nl;
}


void testOverlaps
(
    const treeBoundBox& bb,
    const point& sample,
    const scalar nearestDistSqr
)
{
    FixedList<bool, 8> overlaps;

    for (direction octant = 0; octant < 8; ++octant)
    {
        overlaps[octant] = bb.subOverlaps(octant, sample, nearestDistSqr);
    }

    Info<< "box " << bb << " and "
        << sample << " distSqr:" << nearestDistSqr << nl;

    Info<< "overlaps any:" << bb.overlaps(sample, nearestDistSqr)
        << " octants: " << overlaps << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    treeBoundBox bb(zero_one{});
    treeBoundBox sub(cube(0.1, 0.8));

    Info<< nl
        << "box: " << bb << nl;

    Info<< nl;
    for (direction octant = 0; octant < 8; ++octant)
    {
        Info<< "octant:" << octant
            << " (" << octantName(octant) << ") = "
            << bb.subBbox(octant) << nl;
    }

    Info<< nl;
    for (direction facei = 0; facei < 6; ++facei)
    {
        Info<< "sub-half:" << facei
            << " (" << faceName(facei) << ") = "
            << bb.subHalf(facei) << nl;
    }

    Info<< nl;
    for (direction octant = 0; octant < 8; ++octant)
    {
        const point pt = sub.corner(octant);
        const direction subOctant = bb.subOctant(pt);

        Info<< "point:" << pt
            << " in octant " << subOctant
            << " sub-box: " << bb.subBbox(subOctant) << nl;
    }

    for (const scalar dist : {0.1})
    {
        Info<< nl;
        for (direction octant = 0; octant < 8; ++octant)
        {
            treeBoundBox searchBox(cubeAt(bb.corner(octant), dist));
            testOverlaps(bb, searchBox);
        }

        Info<< nl;
        for (direction facei = 0; facei < 6; ++facei)
        {
            treeBoundBox searchBox(cubeAt(bb.faceCentre(facei), dist));
            testOverlaps(bb, searchBox);
        }
    }

    {
        treeBoundBox largerBox(bb);
        largerBox.grow(0.2);

        // Checking at corners,
        // larger by 0.2 in three directions: radius = 0.3464

        for (const scalar dist : {0.1, 0.35})
        {
            const scalar distSqr = sqr(dist);

            Info<< nl;
            for (direction octant = 0; octant < 8; ++octant)
            {
                testOverlaps(bb, largerBox.corner(octant), distSqr);
            }
        }

        // Checking at face centres,
        // larger by 0.2 in a single direction

        for (const scalar dist : {0.1, 0.25})
        {
            const scalar distSqr = sqr(dist);

            Info<< nl;
            for (direction facei = 0; facei < 6; ++facei)
            {
                testOverlaps(bb, largerBox.faceCentre(facei), distSqr);
            }
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
