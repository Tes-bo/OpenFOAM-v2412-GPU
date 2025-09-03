/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "cell.H"
#include "pyramid.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::cell::typeName = "cell";


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::cell::labels(const faceUList& meshFaces) const
{
    const labelList& cFaces = *this;

    label nVerts = 0;
    for (const label facei : cFaces)
    {
        nVerts += meshFaces[facei].size();
    }

    labelList pointLabels(nVerts);

    // The first face has no duplicates, can copy in values
    const labelList& firstFace = meshFaces[cFaces[0]];

    std::copy(firstFace.cbegin(), firstFace.cend(), pointLabels.begin());

    // Now already contains some vertices
    nVerts = firstFace.size();

    // For the rest of the faces. For each vertex, check if the point is
    // already inserted (up to nVerts, which now carries the number of real
    // points. If not, add it at the end of the list.

    for (label facei = 1; facei < cFaces.size(); ++facei)
    {
        for (const label curPoint : meshFaces[cFaces[facei]])
        {
            bool pointFound = false;

            for (label checki = 0; checki < nVerts; ++checki)
            {
                if (curPoint == pointLabels[checki])
                {
                    pointFound = true;
                    break;
                }
            }

            if (!pointFound)
            {
                pointLabels[nVerts] = curPoint;
                ++nVerts;
            }
        }
    }

    pointLabels.resize(nVerts);

    return pointLabels;
}


Foam::pointField Foam::cell::points
(
    const faceUList& meshFaces,
    const UList<point>& meshPoints
) const
{
    const labelList pointLabels = labels(meshFaces);

    pointField allPoints(pointLabels.size());

    forAll(allPoints, i)
    {
        allPoints[i] = meshPoints[pointLabels[i]];
    }

    return allPoints;
}


Foam::edgeList Foam::cell::edges(const faceUList& meshFaces) const
{
    const labelList& cFaces = *this;

    label nEdges = 0;
    for (const label facei : cFaces)
    {
        nEdges += meshFaces[facei].nEdges();
    }

    edgeList allEdges(nEdges);

    nEdges = 0;

    forAll(cFaces, facei)
    {
        for (const edge& curEdge : meshFaces[cFaces[facei]].edges())
        {
            bool edgeFound = false;

            for (label checki = 0; checki < nEdges; ++checki)
            {
                if (curEdge == allEdges[checki])
                {
                    edgeFound = true;
                    break;
                }
            }

            if (!edgeFound)
            {
                allEdges[nEdges] = curEdge;
                ++nEdges;
            }
        }
    }

    allEdges.resize(nEdges);

    return allEdges;
}


Foam::point Foam::cell::centre
(
    const UList<point>& meshPoints,
    const faceUList& meshFaces
) const
{
    // When one wants to access the cell centre and magnitude, the
    // functionality on the mesh level should be used in preference to the
    // functions provided here. They do not rely to the functionality
    // implemented here, provide additional checking and are more efficient.
    // The cell::centre and cell::mag functions may be removed in the future.

    // WARNING!
    // In the old version of the code, it was possible to establish whether any
    // of the pyramids had a negative volume, caused either by the poor
    // estimate of the cell centre or by the fact that the cell was defined
    // inside out. In the new description of the cell, this can only be
    // established with the reference to the owner-neighbour face-cell
    // relationship, which is not available on this level. Thus, all the
    // pyramids are believed to be positive with no checking.

    // Approximate cell centre as the area average of all face centres

    vector ctr = Zero;
    scalar sumArea = 0;

    const labelList& cFaces = *this;

    for (const label facei : cFaces)
    {
        const scalar magArea = meshFaces[facei].mag(meshPoints);
        ctr += meshFaces[facei].centre(meshPoints)*magArea;
        sumArea += magArea;
    }

    ctr /= sumArea + VSMALL;

    // Calculate the centre by breaking the cell into pyramids and
    // volume-weighted averaging their centres

    scalar sumV = 0;
    vector sumVc = Zero;

    for (const label facei : cFaces)
    {
        const face& f = meshFaces[facei];

        scalar pyrVol = pyramidPointFaceRef(f, ctr).mag(meshPoints);

        // if pyramid inside-out because face points inwards invert
        // N.B. pyramid remains unchanged
        if (pyrVol < 0)
        {
            pyrVol = -pyrVol;
        }

        sumV += pyrVol;
        sumVc += pyrVol * pyramidPointFaceRef(f, ctr).centre(meshPoints);
    }

    return sumVc/(sumV + VSMALL);
}


Foam::scalar Foam::cell::mag
(
    const UList<point>& meshPoints,
    const faceUList& meshFaces
) const
{
    // When one wants to access the cell centre and magnitude, the
    // functionality on the mesh level should be used in preference to the
    // functions provided here. They do not rely to the functionality
    // implemented here, provide additional checking and are more efficient.
    // The cell::centre and cell::mag functions may be removed in the future.

    // WARNING! See cell::centre

    const labelList& cFaces = *this;

    // Approximate cell centre as the average of all face centres

    vector ctr = Zero;
    for (const label facei : cFaces)
    {
        ctr += meshFaces[facei].centre(meshPoints);
    }
    ctr /= cFaces.size();

    // Calculate the magnitude by summing the mags of the pyramids
    scalar sumV = 0;

    for (const label facei : cFaces)
    {
        const face& f = meshFaces[facei];

        sumV += ::Foam::mag(pyramidPointFaceRef(f, ctr).mag(meshPoints));
    }

    return sumV;
}


Foam::Pair<Foam::point>
Foam::cell::box
(
    const UList<point>& meshPoints,
    const faceUList& meshFaces
) const
{
    Pair<point> bb(point::rootMax, point::rootMin);

    for (const label facei : *this)
    {
        for (const label pointi : meshFaces[facei])
        {
            const point& p = meshPoints[pointi];

            bb.first()  = min(bb.first(), p);
            bb.second() = max(bb.second(), p);
        }
    }

    return bb;
}


Foam::Pair<Foam::point> Foam::cell::box(const primitiveMesh& mesh) const
{
    return cell::box(mesh.points(), mesh.faces());
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

bool Foam::operator==(const cell& a, const cell& b)
{
    // Trivial reject: faces are different size
    if (a.size() != b.size())
    {
        return false;
    }

    List<bool> fnd(a.size(), false);

    for (const label curLabel : b)
    {
        bool found = false;

        forAll(a, ai)
        {
            if (a[ai] == curLabel)
            {
                found = true;
                fnd[ai] = true;
                break;
            }
        }

        if (!found)
        {
            return false;
        }
    }

    // Any faces missed?
    forAll(fnd, ai)
    {
        if (!fnd[ai])
        {
            return false;
        }
    }

    return true;
}


// ************************************************************************* //
