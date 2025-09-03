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

#include "collapseEdge.H"

static void markPointNbrs
(
    const triSurface& surf,
    const label facei,
    const bool val,
    boolList& okToCollapse
)
{
    const triSurface::face_type& f = surf.localFaces()[facei];

    forAll(f, fp)
    {
        const labelList& pFaces = surf.pointFaces()[f[fp]];

        forAll(pFaces, i)
        {
            okToCollapse[pFaces[i]] = false;
        }
    }
}


static triSurface pack
(
    const triSurface& surf,
    const pointField& localPoints,
    const labelList& pointMap
)
{
    List<labelledTri> newTriangles(surf.size());
    label nNewTris = 0;

    // Iterate and work on a copy
    for (labelledTri f : surf.localFaces())
    {
        // inplace renumber
        f[0] = pointMap[f[0]];
        f[1] = pointMap[f[1]];
        f[2] = pointMap[f[2]];

        if (f.good())
        {
            newTriangles[nNewTris++] = f;
        }
    }
    newTriangles.resize(nNewTris);

    return triSurface(newTriangles, surf.patches(), localPoints);
}


// Collapses small edge to point, thus removing triangle.
label collapseEdge(triSurface& surf, const scalar minLen)
{
    label nTotalCollapsed = 0;

    while (true)
    {
        const pointField& localPoints = surf.localPoints();
        const List<labelledTri>& localFaces = surf.localFaces();


        // Mapping from old to new points
        labelList pointMap(identity(surf.nPoints()));

        // Storage for new points.
        pointField newPoints(localPoints);

        // To protect neighbours of collapsed faces.
        boolList okToCollapse(surf.size(), true);
        label nCollapsed = 0;

        forAll(localFaces, facei)
        {
            if (okToCollapse[facei])
            {
                // Check edge lengths.
                const triSurface::face_type& f = localFaces[facei];

                forAll(f, fp)
                {
                    label v = f[fp];
                    label v1 = f[f.fcIndex(fp)];

                    if (mag(localPoints[v1] - localPoints[v]) < minLen)
                    {
                        // Collapse f[fp1] onto f[fp].
                        pointMap[v1] = v;
                        newPoints[v] = 0.5*(localPoints[v1] + localPoints[v]);

                        //Pout<< "Collapsing triangle " << facei
                        //    << " to edge mid " << newPoints[v] << endl;

                        nCollapsed++;
                        okToCollapse[facei] = false;

                        // Protect point neighbours from collapsing.
                        markPointNbrs(surf, facei, false, okToCollapse);

                        break;
                    }
                }
            }
        }

        Info<< "collapseEdge : collapsing " << nCollapsed
            << " triangles to a single edge."
            << endl;

        nTotalCollapsed += nCollapsed;

        if (nCollapsed == 0)
        {
            break;
        }

        // Pack the triangles
        surf = pack(surf, newPoints, pointMap);
    }

    // Remove any unused vertices
    surf = triSurface(surf.localFaces(), surf.patches(), surf.localPoints());

    return nTotalCollapsed;
}


// ************************************************************************* //
