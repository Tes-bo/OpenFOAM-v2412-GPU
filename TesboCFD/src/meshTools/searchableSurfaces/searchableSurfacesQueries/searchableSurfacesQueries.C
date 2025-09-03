/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "searchableSurfacesQueries.H"
#include "ListOps.H"
#include "OFstream.H"
#include "meshTools.H"
#include "DynamicField.H"
#include "pointConstraint.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfacesQueries, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfacesQueries::mergeHits
(
    const point& start,

    const label testI,                    // index of surface
    const List<pointIndexHit>& surfHits,  // hits on surface

    labelList& allSurfaces,
    List<pointIndexHit>& allInfo,
    scalarList& allDistSqr
)
{
    // Given current set of hits (allSurfaces, allInfo) merge in those coming
    // from surface surfI.

    // Precalculate distances
    scalarList surfDistSqr(surfHits.size());
    forAll(surfHits, i)
    {
        surfDistSqr[i] = surfHits[i].hitPoint().distSqr(start);
    }

    forAll(surfDistSqr, i)
    {
        label index = findLower(allDistSqr, surfDistSqr[i]);

        // Check if equal to lower.
        if (index >= 0)
        {
            // Same. Do not count.
            //Pout<< "point:" << surfHits[i].point()
            //    << " considered same as:" << allInfo[index].point()
            //    << " within tol:" << mergeDist
            //    << endl;
        }
        else
        {
            // Check if equal to higher
            label next = index + 1;

            if (next < allDistSqr.size())
            {
                //Pout<< "point:" << surfHits[i].point()
                //    << " considered same as:" << allInfo[next].point()
                //    << " within tol:" << mergeDist
                //    << endl;
            }
            else
            {
                // Insert after index
                label sz = allSurfaces.size();
                allSurfaces.setSize(sz+1);
                allInfo.setSize(allSurfaces.size());
                allDistSqr.setSize(allSurfaces.size());
                // Make space.
                for (label j = sz-1; j > index; --j)
                {
                    allSurfaces[j+1] = allSurfaces[j];
                    allInfo[j+1] = allInfo[j];
                    allDistSqr[j+1] = allDistSqr[j];
                }
                // Insert new value
                allSurfaces[index+1] = testI;
                allInfo[index+1] = surfHits[i];
                allDistSqr[index+1] = surfDistSqr[i];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Find any intersection
void Foam::searchableSurfacesQueries::findAnyIntersection
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,
    labelList& hitSurfaces,
    List<pointIndexHit>& hitInfo
)
{
    hitSurfaces.setSize(start.size());
    hitSurfaces = -1;
    hitInfo.setSize(start.size());

    // Work arrays
    labelList hitMap(identity(start.size()));
    pointField p0(start);
    pointField p1(end);
    List<pointIndexHit> intersectInfo(start.size());

    forAll(surfacesToTest, testI)
    {
        // Do synchronised call to all surfaces.
        allSurfaces[surfacesToTest[testI]].findLineAny(p0, p1, intersectInfo);

        // Copy all hits into arguments, continue with misses
        label newI = 0;
        forAll(intersectInfo, i)
        {
            if (intersectInfo[i].hit())
            {
                hitInfo[hitMap[i]] = intersectInfo[i];
                hitSurfaces[hitMap[i]] = testI;
            }
            else
            {
                if (i != newI)
                {
                    hitMap[newI] = hitMap[i];
                    p0[newI] = p0[i];
                    p1[newI] = p1[i];
                }
                newI++;
            }
        }

        // All done? Note that this decision should be synchronised
        if (newI == 0)
        {
            break;
        }

        // Trim and continue
        hitMap.setSize(newI);
        p0.setSize(newI);
        p1.setSize(newI);
        intersectInfo.setSize(newI);
    }
}


void Foam::searchableSurfacesQueries::findAllIntersections
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,
    labelListList& hitSurfaces,
    List<List<pointIndexHit>>& hitInfo
)
{
    // Note: maybe move the single-surface all intersections test into
    // searchable surface? Some of the tolerance issues might be
    // lessened.

    // 2. Currently calling searchableSurface::findLine with start==end
    //    is expected to find no intersection. Problem if it does.

    hitSurfaces.setSize(start.size());
    hitInfo.setSize(start.size());

    if (surfacesToTest.empty())
    {
        return;
    }

    // Test first surface
    allSurfaces[surfacesToTest[0]].findLineAll(start, end, hitInfo);

    // Set hitSurfaces and distance
    List<scalarList> hitDistSqr(hitInfo.size());
    forAll(hitInfo, pointi)
    {
        const List<pointIndexHit>& pHits = hitInfo[pointi];

        labelList& pSurfaces = hitSurfaces[pointi];
        pSurfaces.setSize(pHits.size());
        pSurfaces = 0;

        scalarList& pDistSqr = hitDistSqr[pointi];
        pDistSqr.setSize(pHits.size());
        forAll(pHits, i)
        {
            pDistSqr[i] = pHits[i].hitPoint().distSqr(start[pointi]);
        }
    }


    if (surfacesToTest.size() > 1)
    {
        // Test the other surfaces and merge (according to distance from start).
        for (label testI = 1; testI < surfacesToTest.size(); testI++)
        {
            List<List<pointIndexHit>> surfHits;
            allSurfaces[surfacesToTest[testI]].findLineAll
            (
                start,
                end,
                surfHits
            );

            forAll(surfHits, pointi)
            {
                mergeHits
                (
                    start[pointi],          // Current segment

                    testI,                  // Surface and its hits
                    surfHits[pointi],

                    hitSurfaces[pointi],    // Merge into overall hit info
                    hitInfo[pointi],
                    hitDistSqr[pointi]
                );
            }
        }
    }
}


void Foam::searchableSurfacesQueries::findNearestIntersection
(
   const PtrList<searchableSurface>& allSurfaces,
   const labelList& surfacesToTest,
   const pointField& start,
   const pointField& end,
   labelList& surface1,
   List<pointIndexHit>& hit1,
   labelList& surface2,
   List<pointIndexHit>& hit2
)
{
   // 1. intersection from start to end
   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // Initialize arguments
   surface1.setSize(start.size());
   surface1 = -1;
   hit1.setSize(start.size());

   // Current end of segment to test.
   pointField nearest(end);
   // Work array
   List<pointIndexHit> nearestInfo(start.size());

   forAll(surfacesToTest, testI)
   {
       // See if any intersection between start and current nearest
       allSurfaces[surfacesToTest[testI]].findLine
       (
           start,
           nearest,
           nearestInfo
       );

       forAll(nearestInfo, pointi)
       {
           if (nearestInfo[pointi].hit())
           {
               hit1[pointi] = nearestInfo[pointi];
               surface1[pointi] = testI;
               nearest[pointi] = hit1[pointi].point();
           }
       }
   }


   // 2. intersection from end to last intersection
   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // Find the nearest intersection from end to start. Note that we
   // initialize to the first intersection (if any).
   surface2 = surface1;
   hit2 = hit1;

   // Set current end of segment to test.
   forAll(nearest, pointi)
   {
       if (hit1[pointi].hit())
       {
           nearest[pointi] = hit1[pointi].point();
       }
       else
       {
           // Disable testing by setting to end.
           nearest[pointi] = end[pointi];
       }
   }

   forAll(surfacesToTest, testI)
   {
       // See if any intersection between end and current nearest
       allSurfaces[surfacesToTest[testI]].findLine(end, nearest, nearestInfo);

       forAll(nearestInfo, pointi)
       {
           if (nearestInfo[pointi].hit())
           {
               hit2[pointi] = nearestInfo[pointi];
               surface2[pointi] = testI;
               nearest[pointi] = hit2[pointi].point();
           }
       }
   }
}


void Foam::searchableSurfacesQueries::findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearestSurfaces,
    List<pointIndexHit>& nearestInfo
)
{
    // Find nearest. Return -1 or nearest point

    if (samples.size() != nearestDistSqr.size())
    {
        FatalErrorInFunction << "Inconsistent sizes. samples:" << samples.size()
            << " search-radius:" << nearestDistSqr.size()
            << exit(FatalError);
    }

    // Initialise
    nearestSurfaces.setSize(samples.size());
    nearestSurfaces = -1;
    nearestInfo.setSize(samples.size());

    // Work arrays
    scalarField minDistSqr(nearestDistSqr);
    List<pointIndexHit> hitInfo(samples.size());

    forAll(surfacesToTest, testI)
    {
        allSurfaces[surfacesToTest[testI]].findNearest
        (
            samples,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, pointi)
        {
            if (hitInfo[pointi].hit())
            {
                minDistSqr[pointi] =
                    hitInfo[pointi].point().distSqr(samples[pointi]);

                nearestInfo[pointi] = hitInfo[pointi];
                nearestSurfaces[pointi] = testI;
            }
        }
    }
}


void Foam::searchableSurfacesQueries::findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const labelListList& regionIndices,

    const pointField& samples,
    const scalarField& nearestDistSqr,

    labelList& nearestSurfaces,
    List<pointIndexHit>& nearestInfo
)
{
    // Find nearest. Return -1 or nearest point

    if (samples.size() != nearestDistSqr.size())
    {
        FatalErrorInFunction << "Inconsistent sizes. samples:" << samples.size()
            << " search-radius:" << nearestDistSqr.size()
            << exit(FatalError);
    }


    if (regionIndices.empty())
    {
        findNearest
        (
            allSurfaces,
            surfacesToTest,
            samples,
            nearestDistSqr,
            nearestSurfaces,
            nearestInfo
        );
    }

    // Initialise
    nearestSurfaces.setSize(samples.size());
    nearestSurfaces = -1;
    nearestInfo.setSize(samples.size());

    // Work arrays
    scalarField minDistSqr(nearestDistSqr);
    List<pointIndexHit> hitInfo(samples.size());

    forAll(surfacesToTest, testI)
    {
        allSurfaces[surfacesToTest[testI]].findNearest
        (
            samples,
            minDistSqr,
            regionIndices[testI],
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, pointi)
        {
            if (hitInfo[pointi].hit())
            {
                minDistSqr[pointi] =
                    hitInfo[pointi].point().distSqr(samples[pointi]);

                nearestInfo[pointi] = hitInfo[pointi];
                nearestSurfaces[pointi] = testI;
            }
        }
    }
}


void Foam::searchableSurfacesQueries::findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& start,
    const scalarField& distSqr,
    pointField& near,
    List<pointConstraint>& constraint,
    const label nIter
)
{
    // Multi-surface findNearest


    if (start.size() != distSqr.size())
    {
        FatalErrorInFunction << "Inconsistent sizes. samples:" << start.size()
            << " search-radius:" << distSqr.size()
            << exit(FatalError);
    }


    vectorField normal;
    List<pointIndexHit> info;

    allSurfaces[surfacesToTest[0]].findNearest(start, distSqr, info);
    allSurfaces[surfacesToTest[0]].getNormal(info, normal);

    // Extract useful info from initial start point
    near = start;
    forAll(info, i)
    {
        if (info[i].hit())
        {
            near[i] = info[i].point();
        }
    }

    // Store normal as constraint
    constraint.setSize(near.size());
    constraint = pointConstraint();
    forAll(constraint, i)
    {
        if (info[i].hit())
        {
            constraint[i].applyConstraint(normal[i]);
        }
    }

    if (surfacesToTest.size() >= 2)
    {
        // Work space
        //pointField near1;
        vectorField normal1;

        label surfi = 1;
        for (label iter = 0; iter < nIter; iter++)
        {
            // Find nearest on next surface
            const searchableSurface& s = allSurfaces[surfacesToTest[surfi]];

            // Update: info, normal1
            s.findNearest(near, distSqr, info);
            s.getNormal(info, normal1);

            // Move to intersection of
            //    - previous surface(s) : near+normal
            //    - current surface     : info+normal1
            forAll(near, i)
            {
                if (info[i].hit())
                {
                    if (normal[i] != vector::zero)
                    {
                        // Have previous hit. Find intersection
                        if (mag(normal[i]&normal1[i]) < 1.0-1e-6)
                        {
                            plane pl0(near[i], normal[i], false);
                            plane pl1(info[i].point(), normal1[i], false);

                            plane::ray r(pl0.planeIntersect(pl1));
                            vector n = r.dir() / mag(r.dir());

                            // Calculate vector to move onto intersection line
                            vector d(r.refPoint()-near[i]);
                            d.removeCollinear(n);

                            // Trim the max distance
                            scalar magD = mag(d);
                            if (magD > SMALL)
                            {
                                scalar maxDist = Foam::sqrt(distSqr[i]);
                                if (magD > maxDist)
                                {
                                    // Clip
                                    d /= magD;
                                    d *= maxDist;
                                }

                                near[i] += d;
                                normal[i] = normal1[i];
                                constraint[i].applyConstraint(normal1[i]);
                            }
                        }
                    }
                    else
                    {
                        // First hit
                        near[i] = info[i].point();
                        normal[i] = normal1[i];
                        constraint[i].applyConstraint(normal1[i]);
                    }
                }
            }

            // Step to next surface
            surfi = surfacesToTest.fcIndex(surfi);
        }
    }
}


void Foam::searchableSurfacesQueries::signedDistance
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const volumeType illegalHandling,
    labelList& nearestSurfaces,
    scalarField& distance
)
{
    // Initialise
    distance.setSize(samples.size());
    distance = -GREAT;

    // Find nearest
    List<pointIndexHit> nearestInfo;
    findNearest
    (
        allSurfaces,
        surfacesToTest,
        samples,
        nearestDistSqr,
        nearestSurfaces,
        nearestInfo
    );

    // Determine sign of nearest. Sort by surface to do this.
    DynamicField<point> surfPoints(samples.size());
    DynamicList<label> surfIndices(samples.size());

    forAll(surfacesToTest, testI)
    {
        // Extract samples on this surface
        surfPoints.clear();
        surfIndices.clear();
        forAll(nearestSurfaces, i)
        {
            if (nearestSurfaces[i] == testI)
            {
                surfPoints.append(samples[i]);
                surfIndices.append(i);
            }
        }

        // Calculate sideness of these surface points
        List<volumeType> volType;
        allSurfaces[surfacesToTest[testI]].getVolumeType(surfPoints, volType);

        // Push back to original
        forAll(volType, i)
        {
            label pointi = surfIndices[i];
            scalar dist = samples[pointi].dist(nearestInfo[pointi].hitPoint());

            volumeType vT = volType[i];

            if (vT == volumeType::OUTSIDE)
            {
                distance[pointi] = dist;
            }
            else if (vT == volumeType::INSIDE)
            {
                distance[i] = -dist;
            }
            else
            {
                switch (illegalHandling)
                {
                    case volumeType::OUTSIDE:
                    {
                        distance[pointi] = dist;
                        break;
                    }
                    case volumeType::INSIDE:
                    {
                        distance[pointi] = -dist;
                        break;
                    }
                    default:
                    {
                        FatalErrorInFunction
                            << "getVolumeType failure,"
                            << " neither INSIDE or OUTSIDE."
                            << " point:" << surfPoints[i]
                            << " surface:"
                            << allSurfaces[surfacesToTest[testI]].name()
                            << " volType:" << vT.str()
                            << exit(FatalError);
                        break;
                    }
                }
            }
        }
    }
}


Foam::boundBox Foam::searchableSurfacesQueries::bounds
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelUList& surfacesToTest
)
{
    boundBox bb;

    for (const label surfi : surfacesToTest)
    {
        bb.add(allSurfaces[surfi].bounds());
    }

    return bb;
}


// ************************************************************************* //
