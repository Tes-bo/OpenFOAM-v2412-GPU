/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "cut.H"
#include "addToRunTimeSelectionTable.H"
#include "searchableSurfaces.H"
#include "triSurfaceMesh.H"
#include "searchableBox.H"
#include "searchableRotatedBox.H"
#include "surfaceIntersection.H"
#include "intersectedSurface.H"
#include "edgeIntersections.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace searchableSurfaceModifiers
{
    defineTypeNameAndDebug(cut, 0);
    addToRunTimeSelectionTable(searchableSurfaceModifier, cut, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceModifiers::cut::triangulate
(
    const faceList& fcs,
    pointField& pts,
    triSurface& cutSurf
) const
{
    label nTris = 0;
    forAll(fcs, i)
    {
        nTris += fcs[i].size()-2;
    }

    DynamicList<labelledTri> tris(nTris);

    forAll(fcs, i)
    {
        const face& f = fcs[i];
        // Triangulate around vertex 0
        for (label fp = 1; fp < f.size()-1; fp++)
        {
            tris.append(labelledTri(f[0], f[fp], f[f.fcIndex(fp)], i));
        }
    }
    geometricSurfacePatchList patches(fcs.size());
    forAll(patches, patchi)
    {
        patches[patchi] = geometricSurfacePatch
        (
            geometricSurfacePatch::defaultName(patchi),
            patchi
        );
    }
    cutSurf = triSurface(tris, patches, pts, true);
}


Foam::triSurface& Foam::searchableSurfaceModifiers::cut::triangulate
(
    const searchableSurface& cutter,
    triSurface& cutSurf
) const
{
    if (isA<searchableBox>(cutter))
    {
        const searchableBox& bb = refCast<const searchableBox>(cutter);

        pointField pts(bb.points());
        triangulate(treeBoundBox::faces, pts, cutSurf);

        return cutSurf;
    }
    else if (isA<searchableRotatedBox>(cutter))
    {
        const searchableRotatedBox& bb =
            refCast<const searchableRotatedBox>(cutter);

        pointField pts(bb.points());
        triangulate(treeBoundBox::faces, pts, cutSurf);

        return cutSurf;
    }
    else if (isA<triSurfaceMesh>(cutter))
    {
        return const_cast<triSurfaceMesh&>
        (
            refCast<const triSurfaceMesh>(cutter)
        );
    }
    else
    {
        FatalErrorInFunction
            << "Triangulation only supported for triSurfaceMesh, searchableBox"
            << ", not for surface " << cutter.name()
            << " of type " << cutter.type()
            << exit(FatalError);
        return const_cast<triSurfaceMesh&>
        (
            refCast<const triSurfaceMesh>(cutter)
        );
    }
}


// Keep on shuffling surface points until no more degenerate intersections.
// Moves both surfaces and updates set of edge cuts.
bool Foam::searchableSurfaceModifiers::cut::intersectSurfaces
(
    triSurface& surf1,
    edgeIntersections& edgeCuts1,
    triSurface& surf2,
    edgeIntersections& edgeCuts2
) const
{
    bool hasMoved1 = false;
    bool hasMoved2 = false;

    for (label iter = 0; iter < 10; iter++)
    {
        Info<< "Determining intersections of surf1 edges with surf2"
            << " faces" << endl;

        // Determine surface1 edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters1 = 0;
        {
            triSurfaceSearch querySurf2(surf2);

            scalarField surf1PointTol
            (
                1E-6*edgeIntersections::minEdgeLength(surf1)
            );

            // Determine raw intersections
            edgeCuts1 = edgeIntersections
            (
                surf1,
                querySurf2,
                surf1PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points1(surf1.points());

                nIters1 =
                    edgeCuts1.removeDegenerates
                    (
                        5,              // max iterations
                        surf1,
                        querySurf2,
                        surf1PointTol,
                        points1         // work array
                    );

                if (nIters1 != 0)
                {
                    // Update geometric quantities
                    surf1.movePoints(points1);
                    hasMoved1 = true;
                }
            }
        }

        Info<< "Determining intersections of surf2 edges with surf1"
            << " faces" << endl;

        label nIters2 = 0;
        {
            triSurfaceSearch querySurf1(surf1);

            scalarField surf2PointTol
            (
                1E-6*edgeIntersections::minEdgeLength(surf2)
            );

            // Determine raw intersections
            edgeCuts2 = edgeIntersections
            (
                surf2,
                querySurf1,
                surf2PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points2(surf2.points());

                nIters2 =
                    edgeCuts2.removeDegenerates
                    (
                        5,              // max iterations
                        surf2,
                        querySurf1,
                        surf2PointTol,
                        points2         // work array
                    );

                if (nIters2 != 0)
                {
                    // Update geometric quantities
                    surf2.movePoints(points2);
                    hasMoved2 = true;
                }
            }
        }


        if (nIters1 == 0 && nIters2 == 0)
        {
            //Info<< "** Resolved all intersections to be proper"
            //    << "edge-face pierce" << endl;
            break;
        }
    }

    //if (hasMoved1)
    //{
    //    fileName newFile("surf1.ftr");
    //    Info<< "Surface 1 has been moved. Writing to " << newFile
    //        << endl;
    //    surf1.write(newFile);
    //}
    //
    //if (hasMoved2)
    //{
    //    fileName newFile("surf2.ftr");
    //    Info<< "Surface 2 has been moved. Writing to " << newFile
    //        << endl;
    //    surf2.write(newFile);
    //}

    return hasMoved1 || hasMoved2;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceModifiers::cut::cut
(
    const searchableSurfaces& geometry,
    const dictionary& dict
)
:
    searchableSurfaceModifier(geometry, dict),
    cutterNames_(dict_.get<wordRes>("cutters"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::searchableSurfaceModifiers::cut::modify
(
    const labelList& regions,
    searchableSurface& geom
) const
{
    triSurface& surf = refCast<triSurfaceMesh>(geom);

    bool changed = false;

    // Find the surfaces to cut with
    for (const wordRe& cutterName : cutterNames_)
    {
        labelList geomIDs = wordRes::matching(cutterName, geometry_.names());

        for (const label geomI : geomIDs)
        {
            const searchableSurface& cutter = geometry_[geomI];

            // Triangulate
            triSurface work;
            triSurface& cutSurf = triangulate(cutter, work);

            // Determine intersection (with perturbation)
            edgeIntersections edge1Cuts;
            edgeIntersections edge2Cuts;
            intersectSurfaces
            (
                surf,
                edge1Cuts,
                cutSurf,
                edge2Cuts
            );


            // Determine intersection edges
            surfaceIntersection inter(surf, edge1Cuts, cutSurf, edge2Cuts);


            // Use intersection edges to cut up faces. (does all the hard work)
            intersectedSurface surf3(surf, true, inter);


            // Mark triangles based on whether they are inside or outside
            List<volumeType> volTypes;
            cutter.getVolumeType(surf3.faceCentres(), volTypes);

            label nInside = 0;
            forAll(volTypes, i)
            {
                if (volTypes[i] == volumeType::INSIDE)
                {
                    ++nInside;
                }
            }

            // Add a patch for inside the box
            if (nInside && surf3.patches().size() > 0)
            {
                geometricSurfacePatchList newPatches(surf3.patches());
                label sz = newPatches.size();
                newPatches.setSize(sz+1);
                newPatches[sz] = geometricSurfacePatch
                (
                    newPatches[sz-1].name() + "_inside",
                    newPatches[sz-1].index(),
                    newPatches[sz-1].geometricType()
                );

                Info<< "Moving " << nInside << " out of " << surf3.size()
                    << " triangles to region "
                    << newPatches[sz].name() << endl;


                List<labelledTri> newTris(surf3);
                forAll(volTypes, i)
                {
                    if (volTypes[i] == volumeType::INSIDE)
                    {
                        newTris[i].region() = sz;
                    }
                }
                pointField newPoints(surf3.points());
                surf = triSurface(newTris, newPatches, newPoints, true);

                changed = true;
            }
        }
    }

    return changed;
}


// ************************************************************************* //
