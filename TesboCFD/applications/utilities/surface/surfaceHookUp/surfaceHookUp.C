/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

Application
    surfaceHookUp

Group
    grpSurfaceUtilities

Description
    Find close open edges and stitches the surface along them

Usage
    - surfaceHookUp hookDistance [OPTION]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "triSurfaceMesh.H"
#include "indexedOctree.H"
#include "treeBoundBox.H"
#include "bitSet.H"
#include "unitConversion.H"
#include "searchableSurfaces.H"
#include "IOdictionary.H"

using namespace Foam;

// Split facei along edgeI at position newPointi
void greenRefine
(
    const triSurface& surf,
    const label facei,
    const label edgeI,
    const label newPointi,
    DynamicList<labelledTri>& newFaces
)
{
    const labelledTri& f = surf.localFaces()[facei];
    const edge& e = surf.edges()[edgeI];

    // Find index of edge in face.

    label fp0 = f.find(e[0]);
    label fp1 = f.fcIndex(fp0);
    label fp2 = f.fcIndex(fp1);

    if (f[fp1] == e[1])
    {
        // Edge oriented like face
        newFaces.append
        (
            labelledTri
            (
                f[fp0],
                newPointi,
                f[fp2],
                f.region()
            )
        );
        newFaces.append
        (
            labelledTri
            (
                newPointi,
                f[fp1],
                f[fp2],
                f.region()
            )
        );
    }
    else
    {
        newFaces.append
        (
            labelledTri
            (
                f[fp2],
                newPointi,
                f[fp1],
                f.region()
            )
        );
        newFaces.append
        (
            labelledTri
            (
                newPointi,
                f[fp0],
                f[fp1],
                f.region()
            )
        );
    }
}


//scalar checkEdgeAngle
//(
//    const triSurface& surf,
//    const label edgeIndex,
//    const label pointIndex,
//    const scalar& angle
//)
//{
//    const edge& e = surf.edges()[edgeIndex];

//    const vector eVec = e.unitVec(surf.localPoints());

//    const labelList& pEdges = surf.pointEdges()[pointIndex];
//
//    forAll(pEdges, eI)
//    {
//        const edge& nearE = surf.edges()[pEdges[eI]];

//        const vector nearEVec = nearE.unitVec(surf.localPoints());

//        const scalar dot = eVec & nearEVec;
//        const scalar minCos = degToRad(angle);

//        if (mag(dot) > minCos)
//        {
//            return false;
//        }
//    }

//    return true;
//}


void createBoundaryEdgeTrees
(
    const PtrList<triSurfaceMesh>& surfs,
    PtrList<indexedOctree<treeDataEdge>>& bEdgeTrees,
    labelListList& treeBoundaryEdges
)
{
    forAll(surfs, surfI)
    {
        const triSurface& surf = surfs[surfI];

        // Boundary edges
        treeBoundaryEdges[surfI] =
            identity(surf.nBoundaryEdges(), surf.nInternalEdges());

        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(surf.localPoints()).extend(rndGen, 1e-4, ROOTVSMALL)
        );

        bEdgeTrees.set
        (
            surfI,
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    surf.edges(),
                    surf.localPoints(),
                    treeBoundaryEdges[surfI]    // selected edges
                ),
                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
       );
    }
}


class findNearestOpSubset
{
    const indexedOctree<treeDataEdge>& tree_;

    DynamicList<label>& shapeMask_;

public:

    findNearestOpSubset
    (
        const indexedOctree<treeDataEdge>& tree,
        DynamicList<label>& shapeMask
    )
    :
        tree_(tree),
        shapeMask_(shapeMask)
    {}

    void operator()
    (
        const labelUList& indices,
        const point& sample,

        scalar& nearestDistSqr,
        label& minIndex,
        point& nearestPoint
    ) const
    {
        const treeDataEdge& shape = tree_.shapes();

        for (const label index : indices)
        {
            const label edgeIndex = shape.objectIndex(index);

            if (shapeMask_.found(edgeIndex))
            {
                continue;
            }

            pointHit nearHit = shape.line(index).nearestDist(sample);

            // Only register hit if closest point is not an edge point
            if (nearHit.hit())
            {
                const scalar distSqr = sqr(nearHit.distance());

                if (distSqr < nearestDistSqr)
                {
                    nearestDistSqr = distSqr;
                    minIndex = index;
                    nearestPoint = nearHit.point();
                }
            }
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Hook surfaces to other surfaces by moving and retriangulating their"
        " boundary edges to match other surface boundary edges"
    );
    argList::noParallel();
    argList::addArgument("hookTolerance", "The point merge tolerance");
    argList::addOption("dict", "file", "Alternative surfaceHookUpDict");

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("surfaceHookUpDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    const IOdictionary dict(dictIO);

    const scalar dist(args.get<scalar>(1));
    const scalar matchTolerance(max(1e-6*dist, SMALL));
    const label maxIters = 100;

    Info<< "Hooking distance = " << dist << endl;

    searchableSurfaces surfs
    (
        IOobject
        (
            "surfacesToHook",
            runTime.constant(),
            "triSurface",
            runTime
        ),
        dict,
        true            // assume single-region names get surface name
    );

    Info<< nl << "Reading surfaces: " << endl;
    forAll(surfs, surfI)
    {
        Info<< incrIndent;
        Info<< nl << indent << "Surface     = " << surfs.names()[surfI] << endl;

        const wordList& regions = surfs[surfI].regions();
        forAll(regions, surfRegionI)
        {
            Info<< incrIndent;
            Info<< indent << "Regions = " << regions[surfRegionI] << endl;
            Info<< decrIndent;
        }
        Info<< decrIndent;
    }

    PtrList<indexedOctree<treeDataEdge>> bEdgeTrees(surfs.size());
    labelListList treeBoundaryEdges(surfs.size());

    List<DynamicList<labelledTri>> newFaces(surfs.size());
    List<DynamicList<point>> newPoints(surfs.size());
    List<bitSet> visitedFace(surfs.size());

    PtrList<triSurfaceMesh> newSurfaces(surfs.size());
    forAll(surfs, surfI)
    {
        const triSurfaceMesh& surf =
            refCast<const triSurfaceMesh>(surfs[surfI]);

        newSurfaces.set
        (
            surfI,
            new triSurfaceMesh
            (
                IOobject
                (
                    "hookedSurface_" + surfs.names()[surfI],
                    runTime.constant(),
                    "triSurface",
                    runTime
                ),
                surf
            )
        );
    }

    label nChanged = 0;
    label nIters = 1;

    do
    {
        Info<< nl << "Iteration = " << nIters++ << endl;
        nChanged = 0;

        createBoundaryEdgeTrees(newSurfaces, bEdgeTrees, treeBoundaryEdges);

        forAll(newSurfaces, surfI)
        {
            const triSurface& newSurf = newSurfaces[surfI];

            newFaces[surfI] = newSurf.localFaces();
            newPoints[surfI] = newSurf.localPoints();
            visitedFace[surfI] = bitSet(newSurf.size(), false);
        }

        forAll(newSurfaces, surfI)
        {
            const triSurface& surf = newSurfaces[surfI];

            List<pointIndexHit> bPointsTobEdges(surf.boundaryPoints().size());
            labelList bPointsHitTree(surf.boundaryPoints().size(), -1);

            const labelListList& pointEdges = surf.pointEdges();

            forAll(bPointsTobEdges, bPointi)
            {
                pointIndexHit& nearestHit = bPointsTobEdges[bPointi];

                const label pointi = surf.boundaryPoints()[bPointi];
                const point& samplePt = surf.localPoints()[pointi];

                const labelList& pEdges = pointEdges[pointi];

                // Add edges connected to the edge to the shapeMask
                DynamicList<label> shapeMask;
                shapeMask.append(pEdges);

                forAll(bEdgeTrees, treeI)
                {
                    const indexedOctree<treeDataEdge>& bEdgeTree =
                        bEdgeTrees[treeI];

                    pointIndexHit currentHit =
                        bEdgeTree.findNearest
                        (
                            samplePt,
                            sqr(dist),
                            findNearestOpSubset
                            (
                                bEdgeTree,
                                shapeMask
                            )
                        );

                    if
                    (
                        currentHit.hit()
                     &&
                        (
                            !nearestHit.hit()
                         ||
                            (
                                currentHit.point().distSqr(samplePt)
                              < nearestHit.point().distSqr(samplePt)
                            )
                        )
                    )
                    {
                        nearestHit = currentHit;
                        bPointsHitTree[bPointi] = treeI;
                    }
                }

                if (nearestHit.hit())
                {
    //                bool rejectEdge =
    //                    checkEdgeAngle
    //                    (
    //                        surf,
    //                        nearestHit.index(),
    //                        pointi,
    //                        30
    //                    );

                    scalar distSqr = nearestHit.point().distSqr(samplePt);

                    if (distSqr > Foam::sqr(dist))
                    {
                        nearestHit.setMiss();
                    }
                }
            }

            forAll(bPointsTobEdges, bPointi)
            {
                const pointIndexHit& eHit = bPointsTobEdges[bPointi];

                if (eHit.hit())
                {
                    const label hitSurfI = bPointsHitTree[bPointi];
                    const triSurface& hitSurf = newSurfaces[hitSurfI];

                    const label eIndex =
                        treeBoundaryEdges[hitSurfI][eHit.index()];
                    const edge& e = hitSurf.edges()[eIndex];

                    const label pointi = surf.boundaryPoints()[bPointi];

                    const labelList& eFaces = hitSurf.edgeFaces()[eIndex];

                    if (eFaces.size() != 1)
                    {
                        WarningInFunction
                            << "Edge is attached to " << eFaces.size()
                            << " faces." << endl;

                        continue;
                    }

                    const label facei = eFaces[0];

                    if (visitedFace[hitSurfI][facei])
                    {
                        continue;
                    }

                    DynamicList<labelledTri> newFacesFromSplit(2);

                    const point& pt = surf.localPoints()[pointi];

                    if
                    (
                        (
                            pt.distSqr(hitSurf.localPoints()[e.start()])
                          < matchTolerance
                        )
                     || (
                            pt.distSqr(hitSurf.localPoints()[e.end()])
                          < matchTolerance
                        )
                    )
                    {
                        continue;
                    }

                    nChanged++;

                    label newPointi = -1;

                    // Keep the points in the same place and move the edge
                    if (hitSurfI == surfI)
                    {
                        newPointi = pointi;
                    }
                    else
                    {
                        newPoints[hitSurfI].append(newPoints[surfI][pointi]);
                        newPointi = newPoints[hitSurfI].size() - 1;
                    }

                    // Split the other face.
                    greenRefine
                    (
                        hitSurf,
                        facei,
                        eIndex,
                        newPointi,
                        newFacesFromSplit
                    );

                    visitedFace[hitSurfI].set(facei);

                    forAll(newFacesFromSplit, newFacei)
                    {
                        const labelledTri& fN = newFacesFromSplit[newFacei];

                        if (newFacei == 0)
                        {
                            newFaces[hitSurfI][facei] = fN;
                        }
                        else
                        {
                            newFaces[hitSurfI].append(fN);
                        }
                    }
                }
            }
        }

        Info<< "    Number of edges split = " << nChanged << endl;

        forAll(newSurfaces, surfI)
        {
            newSurfaces.set
            (
                surfI,
                new triSurfaceMesh
                (
                    IOobject
                    (
                        "hookedSurface_" + surfs.names()[surfI],
                        runTime.constant(),
                        "triSurface",
                        runTime
                    ),
                    triSurface
                    (
                        newFaces[surfI],
                        newSurfaces[surfI].patches(),
                        pointField(newPoints[surfI])
                    )
                )
            );
        }

    } while (nChanged > 0 && nIters <= maxIters);

    Info<< endl;

    forAll(newSurfaces, surfI)
    {
        const triSurfaceMesh& newSurf = newSurfaces[surfI];

        Info<< "Writing hooked surface " << newSurf.searchableSurface::name()
            << endl;

        newSurf.searchableSurface::write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
