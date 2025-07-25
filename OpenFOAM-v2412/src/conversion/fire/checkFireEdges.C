/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "checkFireEdges.H"
#include "polyMesh.H"
#include "edgeHashes.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //! \cond fileScope
    //- Print face information for debugging purposes
    static inline void printFace
    (
        const face& f,
        label faceI,
        const edge& currEdge
    )
    {
        Info<< "face " << faceI << ':';
        forAll(f, fpI)
        {
            Info<< ' ';
            if (f[fpI] == currEdge[0] || f[fpI] == currEdge[1])
            {
                Info<< '_';   // highlight the node
            }
            Info<< f[fpI];
        }
        Info<< endl;
    }
    //! \endcond
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::checkFireEdges
(
    const faceList& faces,
    const labelListList& pointFaces,
    const UList<point>& points
)
{
    label nFailedEdges = 0;
    const bool fullCheck = true;

    Info<< "Checking edges according to AVL/FIRE on-the-fly methodology..."
        << endl;

    labelHashSet strayPoints(100);
    edgeHashSet  failedEdges(100);

    forAll(faces, faceI)
    {
        const face& faceA = faces[faceI];

        forAll(faceA, edgei)
        {
            const edge currEdge = faceA.edge(edgei);

            // all faces attached to the first point
            const labelList& otherFaceIds = pointFaces[currEdge[0]];

            forAll(otherFaceIds, otherI)
            {
                const int otherFaceI = otherFaceIds[otherI];
                const face& faceB = faces[otherFaceI];

                // only check once
                if (otherFaceI <= faceI && !fullCheck)
                {
                    continue;
                }

                // get local edges on the second face
                int other_p0 = -1;
                int other_p1 = -1;
                int size_m1  = faceB.size() - 1;

                forAll(faceB, ptI)
                {
                    if (faceB[ptI] == currEdge[0])
                    {
                        other_p0 = ptI;
                    }

                    if (faceB[ptI] == currEdge[1])
                    {
                        other_p1 = ptI;
                    }
                }

                if
                (
                    // did not have both points - can skip
                    other_p0 == -1
                 || other_p1 == -1
                    // a normal edge
                 || abs(other_p0 - other_p1) == 1
                    // handle wrapping
                 || (other_p0 == 0 && other_p1 == size_m1)
                 || (other_p1 == 0 && other_p0 == size_m1)
                )
                {
                    continue;
                }

                // find the "stray" point
                int stray = -1;
                if (abs(other_p0 - other_p1) == 2)
                {
                    // a normal case
                    stray = (other_p0 + other_p1) / 2;
                }
                else if
                (
                    (other_p0 == 0 && other_p1+1 == size_m1)
                 || (other_p1 == 0 && other_p0+1 == size_m1)
                )
                {
                    stray = size_m1;
                }

                if (stray > 0)
                {
                    strayPoints.set(faceB[stray]);
                }

                failedEdges.set(currEdge);

                ++nFailedEdges;

                Info<< nl
                    << "Broken edge calculated between points  "
                    << currEdge[0] << "  " << currEdge[1] << endl;

                printFace(faceA, faceI, currEdge);
                printFace(faceB, otherFaceI, currEdge);
            }
        }
    }

    if (nFailedEdges)
    {
        Info<< endl;
    }
    Info<< "detected " << nFailedEdges << " edge failures";

    // reduce to the actual number of edges
    nFailedEdges = failedEdges.size();

    // report the locations
    if (nFailedEdges)
    {
        Info<< " over " << nFailedEdges << " edges" << endl;

        Info<< nl
            << "edge points" << nl
            << "~~~~~~~~~~~" << endl;

        for (edge thisEdge : failedEdges)  // Use copy of edge
        {
            if (thisEdge.start() > thisEdge.end())
            {
                thisEdge.flip();
            }

            if (notNull(points))
            {
                forAll(thisEdge, keyI)
                {
                    const label ptI = thisEdge[keyI];
                    Info<< "point " << ptI << ": " << points[ptI] << endl;
                }
            }
            else
            {
                forAll(thisEdge, keyI)
                {
                    const label ptI = thisEdge[keyI];
                    Info<< "point " << ptI << endl;
                }
            }
            Info<< endl;
        }

        Info<< nl
            << "stray points" << nl
            << "~~~~~~~~~~~~" << endl;

        {
            labelList keys = strayPoints.sortedToc();

            if (notNull(points))
            {
                forAll(keys, keyI)
                {
                    const label ptI = keys[keyI];
                    Info<< "stray " << ptI << ": " << points[ptI] << endl;
                }
            }
            else
            {
                forAll(keys, keyI)
                {
                    const label ptI = keys[keyI];
                    Info<< "stray " << ptI << endl;
                }
            }

        }
        Info<< endl;
    }
    else
    {
        Info<< endl;
    }

    return nFailedEdges;
}


Foam::label Foam::checkFireEdges
(
    const faceList& faces,
    const UList<point>& points
)
{
    label nPoints = -1;

    if (notNull(points))
    {
        nPoints = points.size();
    }
    else
    {
        // get the max point addressed
        for (const face& f : faces)
        {
            forAll(f, fp)
            {
                if (nPoints < f[fp])
                {
                    nPoints = f[fp];
                }
            }
        }

        ++nPoints;
    }

    labelListList pointFaces(nPoints);
    invertManyToMany(nPoints, faces, pointFaces);

    return checkFireEdges(faces, pointFaces, points);
}


Foam::label Foam::checkFireEdges(const polyMesh& mesh)
{
    return checkFireEdges(mesh.faces(), mesh.pointFaces(), mesh.points());
}


// ************************************************************************* //
