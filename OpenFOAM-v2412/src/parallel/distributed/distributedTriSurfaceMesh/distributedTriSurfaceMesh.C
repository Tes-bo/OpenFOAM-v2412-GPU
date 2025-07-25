/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "distributedTriSurfaceMesh.H"
#include "mapDistribute.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "triangle.H"
#include "matchPoints.H"
#include "globalIndex.H"
#include "Time.H"

#include "IFstream.H"
#include "decompositionMethod.H"
#include "geomDecomp.H"
#include "vectorList.H"
#include "bitSet.H"
#include "PatchTools.H"
#include "OBJstream.H"
#include "labelBits.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedTriSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        distributedTriSurfaceMesh,
        dict
    );


    //- Combine operator for volume types
    class volumeCombineOp
    {
        public:
        void operator()(volumeType& x, const volumeType& y) const
        {
            if (x == volumeType::MIXED || y == volumeType::MIXED)
            {
                FatalErrorInFunction << "Illegal volume type "
                    << volumeType::names[x]
                    << "," << volumeType::names[y] << exit(FatalError);
            }
            else
            {
                switch (x)
                {
                    case volumeType::UNKNOWN:
                    {
                        if (y == volumeType::INSIDE || y == volumeType::OUTSIDE)
                        {
                            x = y;
                        }
                    }
                    break;
                    case volumeType::INSIDE:
                    {
                        if (y == volumeType::OUTSIDE)
                        {
                            FatalErrorInFunction << "Conflicting volume types "
                                << volumeType::names[x] << ","
                                << volumeType::names[y] << exit(FatalError);
                        }
                    }
                    break;
                    case volumeType::OUTSIDE:
                    {
                        if (y == volumeType::INSIDE)
                        {
                            FatalErrorInFunction << "Conflicting volume types "
                                << volumeType::names[x] << ","
                                << volumeType::names[y] << exit(FatalError);
                        }
                    }
                    break;
                    case volumeType::MIXED:
                    break;
                }
            }
        }
    };


    //typedef Tuple2<Pair<point>, volumeType> NearType;
    //
    ////- Combine operator for volume types
    //struct NearTypeCombineOp
    //{
    //    //const auto& this_;
    //
    //    void operator()(NearType& nearX, const NearType& nearY) const
    //    {
    //        const volumeType& x = nearX.second();
    //        const volumeType& y = nearY.second();
    //
    //        if (x == volumeType::MIXED || y == volumeType::MIXED)
    //        {
    //            FatalErrorInFunction << "Illegal volume type "
    //                << volumeType::names[x]
    //                << "," << volumeType::names[y]
    //                << " nearX:" << nearX << " nearY:" << nearY
    //                << exit(FatalError);
    //        }
    //        else
    //        {
    //            switch (x)
    //            {
    //                case volumeType::UNKNOWN:
    //                {
    //                    if (y == volumeType::INSIDE
    //                    || y == volumeType::OUTSIDE)
    //                    {
    //                        nearX = nearY;
    //                    }
    //                }
    //                break;
    //                case volumeType::INSIDE:
    //                {
    //                    if (y == volumeType::OUTSIDE)
    //                    {
    //                        FatalErrorInFunction
    //                            << "Conflicting volume types "
    //                            << volumeType::names[x] << ","
    //                            << volumeType::names[y]
    //                            << " nearX:" << nearX << " nearY:" << nearY
    //                            << exit(FatalError);
    //                    }
    //                }
    //                break;
    //                case volumeType::OUTSIDE:
    //                {
    //                    if (y == volumeType::INSIDE)
    //                    {
    //                        FatalErrorInFunction
    //                            << "Conflicting volume types "
    //                            << volumeType::names[x] << ","
    //                            << volumeType::names[y]
    //                            << " nearX:" << nearX << " nearY:" << nearY
    //                            << exit(FatalError);
    //                    }
    //                }
    //                break;
    //                case volumeType::MIXED:
    //                break;
    //            }
    //        }
    //    }
    //};


    //- Combine operator for nearest
    typedef Tuple2<pointIndexHit, scalar> nearestAndDist;

    const nearestAndDist nearestZero
    (
        nearestAndDist
        (
            pointIndexHit(),
            -GREAT
        )
    );
    class nearestEqOp
    {
        public:
        void operator()(nearestAndDist& x, const nearestAndDist& y) const
        {
            if (x.first().hit())
            {
                if (y.first().hit() && y.second() < x.second())
                {
                    x = y;
                }
            }
            else if (y.first().hit())
            {
                x = y;
            }
        }
    };
}


const Foam::Enum
<
    Foam::distributedTriSurfaceMesh::distributionType
>
Foam::distributedTriSurfaceMesh::distributionTypeNames_
({
    { distributionType::FOLLOW, "follow" },
    { distributionType::INDEPENDENT, "independent" },
    { distributionType::DISTRIBUTED, "distributed" },
    { distributionType::FROZEN, "frozen" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::distributedTriSurfaceMesh::findLocalInstance
(
    const IOobject& io
)
{
    // Modified findInstance which also looks in parent directory

    word instance
    (
        io.time().findInstance
        (
            io.local(),
            word::null,
            IOobject::READ_IF_PRESENT,
            word::null,             // No stop instance
            !UPstream::parRun()     // Fallback to "constant" for non-parallel
        )
    );

    if (instance.size())
    {
        return instance;
    }


    // The rest of this code is only when parRun == true ...

    // Replicate findInstance operation but now on parent directory

    // Search in parent directory
    fileName parentDir =
    (
        io.rootPath()/io.globalCaseName()
       /io.instance()/io.db().dbDir()/io.local()/io.name()
    );

    if (fileHandler().isDir(parentDir))
    {
        return io.instance();
    }


    const scalar startValue = io.time().timeOutputValue();

    instantList ts = io.time().times();

    label instIndex = ts.size()-1;

    // Backward search for first time that is <= startValue
    for (; instIndex >= 0; --instIndex)
    {
        if (ts[instIndex].value() <= startValue)
        {
            break;
        }
    }

    // Continue searching from here
    for (; instIndex >= 0; --instIndex)
    {
        // Shortcut: if actual directory is the timeName we've already tested it
        if (ts[instIndex].name() == io.instance())
        {
            continue;
        }

        parentDir =
        (
            io.rootPath()/io.globalCaseName()
           /ts[instIndex].name()/io.db().dbDir()/io.local()/io.name()
        );

        if (fileHandler().isDir(parentDir))
        {
            return ts[instIndex].name();
        }
    }

    // times() usually already includes the constant() so would have been
    // checked above. Re-test if
    // - times() is empty. Sometimes this can happen (e.g. decomposePar with
    //   collated)
    // - times()[0] is not constant
    if (!ts.size() || ts[0].name() != io.time().constant())
    {
        // Note. This needs to be a hard-coded constant, rather than the
        // constant function of the time, because the latter points to
        // the case constant directory in parallel cases

        parentDir =
        (
            io.rootPath()/io.globalCaseName()
           /io.time().constant()/io.db().dbDir()/io.local()/io.name()
        );

        if (fileHandler().isDir(parentDir))
        {
            return io.time().constant();
        }
    }

    // FatalErrorInFunction
    //     << "Cannot find directory " << io.local() << " in times " << ts
    //     << exit(FatalError);
    //
    // return word::null;

    // Nothing found in parent?
    // - treat like regular findInstance with "constant" fallback
    return io.time().constant();
}


// Read my additional data from the dictionary
bool Foam::distributedTriSurfaceMesh::read()
{
    // Get bb of all domains.
    procBb_.resize_nocopy(Pstream::nProcs());

    if (dict_.empty())
    {
        // Did not find the distributed version; assume master has loaded the
        // triSurfaceMesh version. Make up some settings.

        const boundBox& localBb = triSurfaceMesh::bounds();

        procBb_[Pstream::myProcNo()] =
            treeBoundBoxList(1, treeBoundBox(localBb));

        dict_.add("bounds", procBb_[Pstream::myProcNo()]);

        // Wanted distribution type
        distType_ = DISTRIBUTED;    //INDEPENDENT;
        dict_.add("distributionType", distributionTypeNames_[distType_]);

        // Merge distance
        mergeDist_ = SMALL;
        dict_.add("mergeDistance", mergeDist_);

        // Force underlying triSurfaceMesh to calculate volume type
        // (is topological walk; does not construct tree)
        surfaceClosed_ = triSurfaceMesh::hasVolumeType();
        Pstream::broadcast(surfaceClosed_);
        dict_.add("closed", surfaceClosed_);

        // Delay calculating outside vol type since constructs tree. Is ok
        // after distributing since then local surfaces much smaller
        //outsideVolType_ = volumeType::UNKNOWN;
        //if (surfaceClosed_)
        //{
        //    point outsidePt(localBb.max()+localBb.centre());
        //    List<volumeType> outsideVolTypes;
        //    triSurfaceMesh::getVolumeType
        //    (
        //        pointField(1, outsidePt),
        //        outsideVolTypes
        //    );
        //    outsideVolType_ = outsideVolTypes[0];
        //}
        //dict_.add("outsideVolumeType", volumeType::names[outsideVolType_]);
    }
    else
    {
        dict_.readEntry("bounds", procBb_[Pstream::myProcNo()]);

        // Wanted distribution type
        distType_ = distributionTypeNames_.get("distributionType", dict_);

        // Merge distance
        dict_.readEntry("mergeDistance", mergeDist_);

        // Distribution type
        surfaceClosed_ = dict_.getOrDefault<bool>("closed", false);

        outsideVolType_ = volumeType::names.getOrDefault
        (
            "outsideVolumeType",
            dict_,
            volumeType::UNKNOWN
        );
    }

    Pstream::allGatherList(procBb_);

    return true;
}


// Is segment fully local?
bool Foam::distributedTriSurfaceMesh::isLocal
(
    const List<treeBoundBox>& myBbs,
    const point& start,
    const point& end
)
{
    forAll(myBbs, bbi)
    {
        if (myBbs[bbi].contains(start) && myBbs[bbi].contains(end))
        {
            return true;
        }
    }
    return false;
}


//void Foam::distributedTriSurfaceMesh::splitSegment
//(
//    const label segmenti,
//    const point& start,
//    const point& end,
//    const treeBoundBox& bb,
//
//    DynamicList<segment>& allSegments,
//    DynamicList<label>& allSegmentMap,
//    DynamicList<label> sendMap
//) const
//{
//    // Work points
//    point clipPt0, clipPt1;
//
//    if (bb.contains(start))
//    {
//        // start within, trim end to bb
//        bool clipped = bb.intersects(end, start, clipPt0);
//
//        if (clipped)
//        {
//            // segment from start to clippedStart passes
//            // through proc.
//            sendMap[proci].append(allSegments.size());
//            allSegmentMap.append(segmenti);
//            allSegments.append(segment(start, clipPt0));
//        }
//    }
//    else if (bb.contains(end))
//    {
//        // end within, trim start to bb
//        bool clipped = bb.intersects(start, end, clipPt0);
//
//        if (clipped)
//        {
//            sendMap[proci].append(allSegments.size());
//            allSegmentMap.append(segmenti);
//            allSegments.append(segment(clipPt0, end));
//        }
//    }
//    else
//    {
//        // trim both
//        bool clippedStart = bb.intersects(start, end, clipPt0);
//
//        if (clippedStart)
//        {
//            bool clippedEnd = bb.intersects(end, clipPt0, clipPt1);
//
//            if (clippedEnd)
//            {
//                // middle part of segment passes through proc.
//                sendMap[proci].append(allSegments.size());
//                allSegmentMap.append(segmenti);
//                allSegments.append(segment(clipPt0, clipPt1));
//            }
//        }
//    }
//}


void Foam::distributedTriSurfaceMesh::distributeSegment
(
    const label segmenti,
    const point& start,
    const point& end,

    DynamicList<segment>& allSegments,
    DynamicList<label>& allSegmentMap,
    List<DynamicList<label>>& sendMap
) const
{
    // 1. Fully local already handled outside. Note: retest is cheap.
    if (isLocal(procBb_[Pstream::myProcNo()], start, end))
    {
        return;
    }


    // 2. If fully inside one other processor, then only need to send
    // to that one processor even if it intersects another. Rare occurrence
    // but cheap to test.
    forAll(procBb_, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const List<treeBoundBox>& bbs = procBb_[proci];

            if (isLocal(bbs, start, end))
            {
                sendMap[proci].append(allSegments.size());
                allSegmentMap.append(segmenti);
                allSegments.append(segment(start, end));
                return;
            }
        }
    }


    // 3. If not contained in single processor send to all intersecting
    // processors.
    forAll(procBb_, proci)
    {
        const List<treeBoundBox>& bbs = procBb_[proci];

        forAll(bbs, bbi)
        {
            const treeBoundBox& bb = bbs[bbi];

            // Scheme a: any processor that intersects the segment gets
            // the whole segment.

            // Intersection point
            point clipPt;

            if (bb.intersects(start, end, clipPt))
            {
                sendMap[proci].append(allSegments.size());
                allSegmentMap.append(segmenti);
                allSegments.append(segment(start, end));
            }

            // Alternative: any processor only gets clipped bit of
            // segment. This gives small problems with additional
            // truncation errors.
            //splitSegment
            //(
            //    segmenti,
            //    start,
            //    end,
            //    bb,
            //
            //    allSegments,
            //    allSegmentMap,
            //   sendMap[proci]
            //);
        }
    }
}


Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::distributeSegments
(
    const pointField& start,
    const pointField& end,

    List<segment>& allSegments,
    labelList& allSegmentMap
) const
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    labelListList sendMap(Pstream::nProcs());

    {
        // Since intersection test is quite expensive compared to memory
        // allocation we use DynamicList to immediately store the segment
        // in the correct bin.

        // Segments to test
        DynamicList<segment> dynAllSegments(start.size());
        // Original index of segment
        DynamicList<label> dynAllSegmentMap(start.size());
        // Per processor indices into allSegments to send
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());

        // Pre-size
        forAll(dynSendMap, proci)
        {
            dynSendMap[proci].reserve
            (
                (proci == Pstream::myProcNo())
              ? start.size()
              : start.size()/Pstream::nProcs()
            );
        }

        forAll(start, segmenti)
        {
            distributeSegment
            (
                segmenti,
                start[segmenti],
                end[segmenti],

                dynAllSegments,
                dynAllSegmentMap,
                dynSendMap
            );
        }

        // Convert dynamicList to labelList
        sendMap.resize_nocopy(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }

        allSegments.transfer(dynAllSegments);
        allSegmentMap.transfer(dynAllSegmentMap);
    }

    return autoPtr<mapDistribute>::New(std::move(sendMap));
}


void Foam::distributedTriSurfaceMesh::findLine
(
    const bool nearestIntersection,
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findLine :"
            << " intersecting surface " << searchableSurface::name()
            << " with "
            << start.size() << " rays" << endl;
    }
    addProfiling(findLine, "distributedTriSurfaceMesh::findLine");

    const indexedOctree<treeDataTriSurface>& octree = tree();

    // Initialise
    info.setSize(start.size());
    forAll(info, i)
    {
        info[i].setMiss();
    }

    // Important:force synchronised construction of indexing
    const globalIndex& triIndexer = globalTris();


    // Do any local queries
    // ~~~~~~~~~~~~~~~~~~~~

    label nLocal = 0;

    forAll(start, i)
    {
        if (isLocal(procBb_[Pstream::myProcNo()], start[i], end[i]))
        {
            if (nearestIntersection)
            {
                info[i] = octree.findLine(start[i], end[i]);
            }
            else
            {
                info[i] = octree.findLineAny(start[i], end[i]);
            }

            if (info[i].hit())
            {
                info[i].setIndex(triIndexer.toGlobal(info[i].index()));
            }
            nLocal++;
        }
    }


    if
    (
        returnReduce(nLocal, sumOp<label>())
      < returnReduce(start.size(), sumOp<label>())
    )
    {
        // Not all can be resolved locally. Build segments and map,
        // send over segments, do intersections, send back and merge.


        // Construct queries (segments)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Segments to test
        List<segment> allSegments(start.size());
        // Original index of segment
        labelList allSegmentMap(start.size());

        const autoPtr<mapDistribute> mapPtr
        (
            distributeSegments
            (
                start,
                end,
                allSegments,
                allSegmentMap
            )
        );
        const mapDistribute& map = mapPtr();

        label nOldAllSegments = allSegments.size();


        // Exchange the segments
        // ~~~~~~~~~~~~~~~~~~~~~

        map.distribute(allSegments);


        // Do tests i need to do
        // ~~~~~~~~~~~~~~~~~~~~~

        // Intersections
        List<pointIndexHit> intersections(allSegments.size());

        forAll(allSegments, i)
        {
            if (nearestIntersection)
            {
                intersections[i] = octree.findLine
                (
                    allSegments[i].first(),
                    allSegments[i].second()
                );
            }
            else
            {
                intersections[i] = octree.findLineAny
                (
                    allSegments[i].first(),
                    allSegments[i].second()
                );
            }

            // Convert triangle index to global numbering
            if (intersections[i].hit())
            {
                intersections[i].setIndex
                (
                    triIndexer.toGlobal(intersections[i].index())
                );
            }
        }


        // Exchange the intersections (opposite to segments)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        map.reverseDistribute(nOldAllSegments, intersections);


        // Extract the hits
        // ~~~~~~~~~~~~~~~~

        forAll(intersections, i)
        {
            const pointIndexHit& allInfo = intersections[i];
            label segmenti = allSegmentMap[i];
            pointIndexHit& hitInfo = info[segmenti];

            if (allInfo.hit())
            {
                if (!hitInfo.hit())
                {
                    // No intersection yet so take this one
                    hitInfo = allInfo;
                }
                else if (nearestIntersection)
                {
                    // Nearest intersection
                    if
                    (
                        start[segmenti].distSqr(allInfo.point())
                      < start[segmenti].distSqr(hitInfo.point())
                    )
                    {
                        hitInfo = allInfo;
                    }
                }
            }
        }
    }
}


void Foam::distributedTriSurfaceMesh::convertTriIndices
(
    List<pointIndexHit>& info
) const
{
    // Important:force synchronised construction of indexing
    const globalIndex& triIndexer = globalTris();

    for (pointIndexHit& pi : info)
    {
        if (pi.hit())
        {
            pi.setIndex(triIndexer.toGlobal(pi.index()));
        }
    }
}


// Exchanges indices to the processor they come from.
// - calculates exchange map
// - uses map to calculate local triangle index
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::calcLocalQueries
(
    const List<pointIndexHit>& info,
    labelList& triangleIndex
) const
{
    // Note: does not filter duplicate queries so a triangle that gets requested
    //       from more than one processor also get local queried more than
    //       once.

    triangleIndex.setSize(info.size());

    const globalIndex& triIndexer = globalTris();


    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // Since determining which processor the query should go to is
    // cheap we do a multi-pass algorithm to save some memory temporarily.

    // 1. Count
    labelList nSend(Pstream::nProcs(), Zero);

    forAll(info, i)
    {
        if (info[i].hit())
        {
            label proci = triIndexer.whichProcID(info[i].index());
            nSend[proci]++;
        }
    }

    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());
    forAll(nSend, proci)
    {
        sendMap[proci].setSize(nSend[proci]);
        nSend[proci] = 0;
    }

    // 3. Fill sendMap
    forAll(info, i)
    {
        if (info[i].hit())
        {
            label proci = triIndexer.whichProcID(info[i].index());
            triangleIndex[i] = triIndexer.toLocal(proci, info[i].index());
            sendMap[proci][nSend[proci]++] = i;
        }
        else
        {
            triangleIndex[i] = -1;
        }
    }


    // Pack into distribution map
    auto mapPtr = autoPtr<mapDistribute>::New(std::move(sendMap));

    // Send over queries
    mapPtr().distribute(triangleIndex);

    return mapPtr;
}


bool Foam::distributedTriSurfaceMesh::contains
(
    const List<treeBoundBox>& bbs,
    const point& sample
) const
{
    forAll(bbs, bbi)
    {
        if (bbs[bbi].contains(sample))
        {
            return true;
        }
    }
    return false;
}


Foam::Tuple2<Foam::label, Foam::scalar>
Foam::distributedTriSurfaceMesh::findBestProcs
(
    const point& centre,
    const scalar radiusSqr,
    boolList& procContains,
    boolList& procOverlaps,
    label& minProci
) const
{
    // Find processors:
    // - that contain the centre
    // - or overlap the search sphere

    procContains.setSize(Pstream::nProcs());
    procContains = false;

    procOverlaps.setSize(Pstream::nProcs());
    procOverlaps = false;

    minProci = -1;

    scalar minDistSqr = radiusSqr;

    label nContain = 0;
    forAll(procBb_, proci)
    {
        const List<treeBoundBox>& bbs = procBb_[proci];

        forAll(bbs, bbi)
        {
            if (bbs[bbi].contains(centre))
            {
                // We found a bb that contains the centre. There must be
                // a triangle inside (since otherwise the bb would never
                // have been created).
                if (!procContains[proci])
                {
                    procContains[proci] = true;
                    nContain++;
                }
                // Minimum search distance to find the triangle
                point near, far;
                bbs[bbi].calcExtremities(centre, near, far);
                minDistSqr = min(minDistSqr, centre.distSqr(far));
            }
        }
    }

    if (nContain > 0)
    {
        return Tuple2<label, scalar>(nContain, minDistSqr);
    }
    else
    {
        scalar maxDistSqr = radiusSqr;

        // Pass 1: find box with nearest minimum distance. Store its maximum
        //         extent as well. Since box will always contain a triangle
        //         this guarantees at least one hit.
        forAll(procBb_, proci)
        {
            const List<treeBoundBox>& bbs = procBb_[proci];

            forAll(bbs, bbi)
            {
                if (bbs[bbi].overlaps(centre, radiusSqr))
                {
                    point near, far;
                    bbs[bbi].calcExtremities(centre, near, far);

                    scalar d2 = centre.distSqr(near);
                    if (d2 < minDistSqr)
                    {
                        minDistSqr = d2;
                        maxDistSqr = min(radiusSqr, centre.distSqr(far));
                        minProci = proci;
                    }
                }
            }
        }

        label nOverlap = 0;
        if (minProci >= 0)
        {
            // Pass 2. Find all bb in range minDistSqr..maxDistSqr

            procOverlaps[minProci] = true;
            nOverlap++;

            forAll(procBb_, proci)
            {
                if (proci != minProci)
                {
                    const List<treeBoundBox>& bbs = procBb_[proci];
                    forAll(bbs, bbi)
                    {
                        if (bbs[bbi].overlaps(centre, maxDistSqr))
                        {
                            procOverlaps[proci] = true;
                            nOverlap++;
                            break;
                        }
                    }
                }
            }
        }
        return Tuple2<label, scalar>(nOverlap, maxDistSqr);
    }
}


Foam::label Foam::distributedTriSurfaceMesh::calcOverlappingProcs
(
    const point& centre,
    const scalar radiusSqr,
    boolList& overlaps
) const
{
    overlaps = false;
    label nOverlaps = 0;

    forAll(procBb_, proci)
    {
        const List<treeBoundBox>& bbs = procBb_[proci];

        forAll(bbs, bbi)
        {
            if (bbs[bbi].overlaps(centre, radiusSqr))
            {
                overlaps[proci] = true;
                nOverlaps++;
                break;
            }
        }
    }
    return nOverlaps;
}


// Generate queries for parallel distance calculation
// - calculates exchange map
// - uses map to exchange points and radius
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::calcLocalQueries
(
    const bool includeLocalProcessor,
    const pointField& centres,
    const scalarField& radiusSqr,

    pointField& allCentres,
    scalarField& allRadiusSqr,
    labelList& allSegmentMap
) const
{
    // Determine queries
    // ~~~~~~~~~~~~~~~~~

    labelListList sendMap(Pstream::nProcs());

    {
        // Queries
        DynamicList<point> dynAllCentres(centres.size());
        DynamicList<scalar> dynAllRadiusSqr(centres.size());
        // Original index of segment
        DynamicList<label> dynAllSegmentMap(centres.size());
        // Per processor indices into allSegments to send
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());

        // Pre-size
        forAll(dynSendMap, proci)
        {
            dynSendMap[proci].reserve
            (
                (proci == Pstream::myProcNo())
              ? centres.size()
              : centres.size()/Pstream::nProcs()
            );
        }

        // Work array - whether processor bb overlaps the bounding sphere.
        boolList procBbOverlaps(Pstream::nProcs());

        forAll(centres, centrei)
        {
            // Find the processor this sample+radius overlaps.
            calcOverlappingProcs
            (
                centres[centrei],
                radiusSqr[centrei],
                procBbOverlaps
            );

            forAll(procBbOverlaps, proci)
            {
                if
                (
                    procBbOverlaps[proci]
                 && (
                        includeLocalProcessor
                     || proci != Pstream::myProcNo()
                    )
                )
                {
                    dynSendMap[proci].append(dynAllCentres.size());
                    dynAllSegmentMap.append(centrei);
                    dynAllCentres.append(centres[centrei]);
                    dynAllRadiusSqr.append(radiusSqr[centrei]);
                }
            }
        }

        // Convert dynamicList to labelList
        sendMap.resize_nocopy(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }

        allCentres.transfer(dynAllCentres);
        allRadiusSqr.transfer(dynAllRadiusSqr);
        allSegmentMap.transfer(dynAllSegmentMap);
    }

    return autoPtr<mapDistribute>::New(std::move(sendMap));
}


Foam::volumeType Foam::distributedTriSurfaceMesh::edgeSide
(
    const point& sample,
    const point& nearestPoint,
    const label face0,
    const label face1
) const
{
    const triSurface& surf = static_cast<const triSurface&>(*this);
    const pointField& points = surf.points();

    // Compare to bisector. This is actually correct since edge is
    // nearest so there is a knife-edge.

    //const vectorField& faceNormals = surf.faceNormals();
    //vector n = faceNormals[face0] + faceNormals[face1];
    vector n = surf[face0].unitNormal(points)+surf[face1].unitNormal(points);

    if (((sample - nearestPoint) & n) > 0)
    {
        return volumeType::OUTSIDE;
    }
    else
    {
        return volumeType::INSIDE;
    }
}


Foam::label Foam::distributedTriSurfaceMesh::findOtherFace
(
    const labelListList& pointFaces,
    const label nearFacei,
    const label nearLabel
) const
{
    const triSurface& surf = static_cast<const triSurface&>(*this);
    const triSurface::face_type& nearF = surf[nearFacei];

    const edge e(nearF[nearLabel], nearF[nearF.fcIndex(nearLabel)]);

    const labelList& pFaces = pointFaces[e[0]];
    for (const label facei : pFaces)
    {
        if (facei != nearFacei)
        {
            if (surf[facei].contains(e))
            {
                return facei;
            }
        }
    }
    return -1;
}


void Foam::distributedTriSurfaceMesh::calcFaceFaces
(
    const triSurface& s,
    const labelListList& pointFaces,
    labelListList& faceFaces
)
{
    faceFaces.setSize(s.size());

    DynamicList<label> nbrs;

    forAll(faceFaces, facei)
    {
        const labelledTri& f = s[facei];

        nbrs.reserve(f.size());
        nbrs.clear();

        forAll(f, fp)
        {
            const edge e(f[fp], f[f.fcIndex(fp)]);
            const labelList& pFaces = pointFaces[f[fp]];
            for (const label otherFacei : pFaces)
            {
                if (otherFacei != facei)
                {
                    if (s[otherFacei].contains(e))
                    {
                        if (!nbrs.find(otherFacei))
                        {
                            nbrs.append(otherFacei);
                        }
                    }
                }
            }
        }
        faceFaces[facei] = std::move(nbrs);
    }
}


void Foam::distributedTriSurfaceMesh::surfaceSide
(
    const pointField& samples,
    const List<pointIndexHit>& nearestInfo,
    List<volumeType>& volType
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::surfaceSide :"
            << " on surface " << searchableSurface::name()
            << " finding surface side given points on surface for "
            << samples.size() << " samples" << endl;
    }

    // Use global index to send local tri and nearest back to originating
    // processor

    labelList triangleIndex(nearestInfo.size());
    autoPtr<mapDistribute> mapPtr
    (
        calcLocalQueries
        (
            nearestInfo,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();

    // Send over samples
    pointField localSamples(samples);
    map.distribute(localSamples);


    // Do my tests
    // ~~~~~~~~~~~

    volType.setSize(triangleIndex.size());
    volType = volumeType::UNKNOWN;

    const triSurface& surf = static_cast<const triSurface&>(*this);
    const pointField& points = surf.points();
    {
        //const labelListList& pointFaces = surf.pointFaces();
        // Construct pointFaces. Let's hope surface has compact point
        // numbering ...
        labelListList pointFaces;
        invertManyToMany(points.size(), surf, pointFaces);

        EdgeMap<labelPair> edgeToFaces;

        forAll(triangleIndex, i)
        {
            const label facei = triangleIndex[i];
            const triSurface::face_type& f = surf[facei];
            const point& sample = localSamples[i];

            // Find where point is on face
            label nearType, nearLabel;
            pointHit pHit =
                f.nearestPointClassify(sample, points, nearType, nearLabel);

            const point& nearestPoint(pHit.point());

            if (nearType == triPointRef::NONE)
            {
                const vector sampleNearestVec = (sample - nearestPoint);

                // Nearest to face interior. Use faceNormal to determine side
                //scalar c = sampleNearestVec & surf.faceNormals()[facei];
                scalar c = sampleNearestVec & surf[facei].areaNormal(points);

                if (c > 0)
                {
                    volType[i] = volumeType::OUTSIDE;
                }
                else
                {
                    volType[i] = volumeType::INSIDE;
                }
            }
            else if (nearType == triPointRef::EDGE)
            {
                // Nearest to edge nearLabel. Note that this can only be a
                // knife-edge
                // situation since otherwise the nearest point could never be
                // the edge.

                label otherFacei = findOtherFace(pointFaces, facei, nearLabel);
                if (otherFacei != -1)
                {
                    volType[i] =
                        edgeSide(sample, nearestPoint, facei, otherFacei);
                }
                else
                {
                    // Open edge. Leave volume type unknown
                }
            }
            else
            {
                // Nearest to point. Could use pointNormal here but is not
                // correct.
                // Instead determine which edge using point is nearest and
                // use test above (nearType == triPointRef::EDGE).

                const label pointi = f[nearLabel];
                const labelList& pFaces = pointFaces[pointi];
                const vector sampleNearestVec = (sample - nearestPoint);

                // Loop over all faces. Check if have both edge faces. Do
                // test.
                edgeToFaces.clear();

                scalar maxCosAngle = -GREAT;
                labelPair maxEdgeFaces(-1, -1);

                for (const label facei : pFaces)
                {
                    const triSurface::face_type& f = surf[facei];

                    label fp = f.find(pointi);
                    label p1 = f[f.fcIndex(fp)];
                    label pMin1 = f[f.rcIndex(fp)];

                    Pair<edge> edges
                    (
                        edge(pointi, p1),
                        edge(pointi, pMin1)
                    );

                    // Check edge fp-to-fp+1  and fp-1
                    // determine distance/angle to nearPoint
                    for (const edge& e : edges)
                    {
                        auto iter = edgeToFaces.find(e);
                        if (iter.good())
                        {
                            if (iter().second() == -1)
                            {
                                // Found second face. Now we have edge+faces
                                iter().second() = facei;

                                vector eVec(e.vec(points));
                                scalar magEVec = mag(eVec);

                                if (magEVec > VSMALL)
                                {
                                    eVec /= magEVec;

                                    // Determine edge most in direction of
                                    // sample
                                    scalar cosAngle(sampleNearestVec&eVec);
                                    if (cosAngle > maxCosAngle)
                                    {
                                        maxCosAngle = cosAngle;
                                        maxEdgeFaces = iter();
                                    }
                                }
                            }
                            else
                            {
                                FatalErrorInFunction << "Not closed"
                                    << exit(FatalError);
                            }
                        }
                        else
                        {
                            edgeToFaces.insert(e, labelPair(facei, -1));
                        }
                    }
                }


                // Check that surface is closed
                bool closed = true;
                for (auto iter : edgeToFaces)
                {
                    if (iter[0] == -1 || iter[1] == -1)
                    {
                        closed = false;
                        break;
                    }
                }
                if (closed)
                {
                    volType[i] = edgeSide
                    (
                        sample,
                        nearestPoint,
                        maxEdgeFaces[0],
                        maxEdgeFaces[1]
                    );
                }
            }
        }
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    // Note that we make sure that if multiple processors hold data only
    // the one with inside/outside wins. (though this should already be
    // handled by the fact we have a unique nearest triangle so we only
    // send the volume-query to a single processor)


    //forAll(localSamples, i)
    //{
    //    Pout<< "surfaceSide : for localSample:" << localSamples[i]
    //        << " found volType:" << volumeType::names[volType[i]]
    //        << endl;
    //}

    const volumeType zero(volumeType::UNKNOWN);
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        nearestInfo.size(),
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        volType,
        zero,
        volumeCombineOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map.comm()
    );

    //// Pack sample, nearest and volType so we get nicer error messages
    ////typedef Tuple2<Pair<point>, volumeType>> NearType;
    //List<NearType> nearTypes(volType.size());
    //forAll(localSamples, i)
    //{
    //    const point& sample = localSamples[i];
    //    const point& near = nearest[i];
    //    nearTypes[i] = NearType(Pair<point>(sample, near), volType[i]);
    //}
    //
    //
    //const NearType zero(Pair<point>(Zero, Zero), volumeType::UNKNOWN);
    //mapDistributeBase::distribute
    //(
    //    Pstream::commsTypes::nonBlocking,
    //    List<labelPair>::null(),
    //    nearestInfo.size(),
    //    map.constructMap(),
    //    map.constructHasFlip(),
    //    map.subMap(),
    //    map.subHasFlip(),
    //    nearTypes,
    //    zero,
    //    NearTypeCombineOp(),
    //    noOp(),           // no flipping
    //    UPstream::msgType(),
    //    map.comm()
    //);
    //volType.setSize(nearTypes.size());
    //forAll(nearTypes, i)
    //{
    //    volType[i] = nearTypes[i].second();
    //}

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::surfaceSide :"
            << " finished finding surface" << searchableSurface::name()
            << " given points on surface "
            << searchableSurface::name() << " for "
            << samples.size() << " samples" << endl;
    }
}


void Foam::distributedTriSurfaceMesh::collectLeafMids
(
    const label nodeI,
    DynamicField<point>& midPoints
) const
{
    const auto& nod = tree().nodes()[nodeI];

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        const labelBits index = nod.subNodes_[octant];

        if (indexedOctree<treeDataTriSurface>::isNode(index))
        {
            // tree node. Recurse.
            collectLeafMids
            (
                indexedOctree<treeDataTriSurface>::getNode(index),
                midPoints
            );
        }
        else if (indexedOctree<treeDataTriSurface>::isContent(index))
        {}
        else
        {
            // No data in this octant. Set type for octant acc. to the mid
            // of its bounding box.
            const treeBoundBox subBb = nod.bb_.subBbox(octant);
            midPoints.append(subBb.centre());
        }
    }
}


Foam::volumeType Foam::distributedTriSurfaceMesh::calcVolumeType
(
    const List<volumeType>& midPointTypes,
    label& midPointi,
    PackedList<2>& nodeTypes,
    const label nodeI
) const
{
    // Pre-calculates wherever possible the volume status per node/subnode.
    // Recurses to determine status of lowest level boxes. Level above is
    // combination of octants below.

    const auto& nod = tree().nodes()[nodeI];

    volumeType myType = volumeType::UNKNOWN;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        volumeType subType;

        const labelBits index = nod.subNodes_[octant];

        if (indexedOctree<treeDataTriSurface>::isNode(index))
        {
            // tree node. Recurse.
            subType = calcVolumeType
            (
                midPointTypes,
                midPointi,
                nodeTypes,
                indexedOctree<treeDataTriSurface>::getNode(index)
            );
        }
        else if (indexedOctree<treeDataTriSurface>::isContent(index))
        {
            // Contents. Depending on position in box might be on either
            // side.
            subType = volumeType::MIXED;
        }
        else
        {
            // No data in this octant. Set type for octant acc. to the mid
            // of its bounding box.
            //Pout<< "    for leaf at bb:" << nod.bb_.subBbox(octant)
            //    << " node:" << nodeI
            //    << " octant:" << octant
            //    << " caching volType to:" << midPointTypes[midPointi] << endl;
            subType = midPointTypes[midPointi++];
        }

        // Store octant type
        nodeTypes.set(labelBits::pack(nodeI, octant), subType);

        // Combine sub node types into type for treeNode. Result is 'mixed' if
        // types differ among subnodes.
        if (myType == volumeType::UNKNOWN)
        {
            myType = subType;
        }
        else if (subType != myType)
        {
            myType = volumeType::MIXED;
        }
    }
    return myType;
}


Foam::volumeType Foam::distributedTriSurfaceMesh::cachedVolumeType
(
    const label nodeI,
    const point& sample
) const
{
    const auto& nod = tree().nodes()[nodeI];

    direction octant = nod.bb_.subOctant(sample);

    volumeType octantType = volumeType::type
    (
        tree().nodeTypes().get(labelBits::pack(nodeI, octant))
    );

    if (octantType == volumeType::INSIDE)
    {
        return octantType;
    }
    else if (octantType == volumeType::OUTSIDE)
    {
        return octantType;
    }
    else if (octantType == volumeType::UNKNOWN)
    {
        // Can happen for e.g. non-manifold surfaces.
        return octantType;
    }
    else if (octantType == volumeType::MIXED)
    {
        labelBits index = nod.subNodes_[octant];

        if (indexedOctree<treeDataTriSurface>::isNode(index))
        {
            // Recurse
            volumeType subType = cachedVolumeType
            (
                indexedOctree<treeDataTriSurface>::getNode(index),
                sample
            );

            return subType;
        }
        else if (indexedOctree<treeDataTriSurface>::isContent(index))
        {
            // Content.
            return volumeType::UNKNOWN;
        }
        else
        {
            // Empty node. Cannot have 'mixed' as its type since not divided
            // up and has no items inside it.
            FatalErrorInFunction
                << "Sample:" << sample << " node:" << nodeI
                << " with bb:" << nod.bb_ << nl
                << "Empty subnode has invalid volume type MIXED."
                << abort(FatalError);

            return volumeType::UNKNOWN;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Sample:" << sample << " at node:" << nodeI
            << " octant:" << octant
            << " with bb:" << nod.bb_.subBbox(octant) << nl
            << "Node has invalid volume type " << octantType
            << abort(FatalError);

        return volumeType::UNKNOWN;
    }
}


// Find bounding boxes that guarantee a more or less uniform distribution
// of triangles. Decomposition in here is only used to get the bounding
// boxes, actual decomposition is done later on.
// Returns a per processor a list of bounding boxes that most accurately
// describe the shape. For now just a single bounding box per processor but
// optimisation might be to determine a better fitting shape.
Foam::List<Foam::List<Foam::treeBoundBox>>
Foam::distributedTriSurfaceMesh::independentlyDistributedBbs
(
    const triSurface& s
)
{
    addProfiling
    (
        distribute,
        "distributedTriSurfaceMesh::independentlyDistributedBbs"
    );

    if (!decomposer_)
    {
        // Use singleton decomposeParDict. Cannot use decompositionModel
        // here since we've only got Time and not a mesh.

        const auto* dictPtr =
            searchableSurface::time().findObject<IOdictionary>
            (
                // == decompositionModel::canonicalName
                "decomposeParDict"
            );

        if (dictPtr)
        {
            decomposer_ = decompositionMethod::New(*dictPtr);
        }
        else
        {
            if (!decomposeParDict_)
            {
                decomposeParDict_.reset
                (
                    new IOdictionary
                    (
                        IOobject
                        (
                            // == decompositionModel::canonicalName
                            "decomposeParDict",
                            searchableSurface::time().system(),
                            searchableSurface::time(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    )
                );
            }
            decomposer_ = decompositionMethod::New(*decomposeParDict_);
        }
    }


    // Initialise to inverted box
    List<List<treeBoundBox>> bbs
    (
        Pstream::nProcs(),
        List<treeBoundBox>(1, treeBoundBox::null())
    );

    const globalIndex& triIndexer = globalTris();

    bool masterOnly;
    {
        masterOnly = true;
        for (const int proci : Pstream::subProcs())
        {
            if (triIndexer.localSize(proci) != 0)
            {
                masterOnly = false;
                break;
            }
        }
    }

    if (masterOnly)
    {
        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
                << " determining master-only decomposition for " << s.size()
                << " centroids for " << searchableSurface::name() << endl;
        }

        // Triangle centres
        pointField triCentres(s.size());
        forAll(s, trii)
        {
            triCentres[trii] = s[trii].centre(s.points());
        }

        labelList distribution;
        if (!isA<geomDecomp>(decomposer_()))
        {
            // Use connectivity
            labelListList pointFaces;
            invertManyToMany(s.points().size(), s, pointFaces);
            labelListList faceFaces(s.size());
            calcFaceFaces(s, pointFaces, faceFaces);

            // Do the actual decomposition
            const bool oldParRun = UPstream::parRun(false);
            distribution = decomposer_().decompose(faceFaces, triCentres);
            UPstream::parRun(oldParRun);  // Restore parallel state
        }
        else
        {
            // Do the actual decomposition
            distribution = decomposer_().decompose(triCentres);
        }

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
                << " determining processor bounding boxes for surface"
                << searchableSurface::name() << endl;
        }

        // Find bounding box for all triangles on new distribution.
        forAll(s, trii)
        {
            treeBoundBox& bb = bbs[distribution[trii]][0];
            bb.add(s.points(), s[trii]);
        }

        // Now combine for all processors and convert to correct format.
        forAll(bbs, proci)
        {
            // Inflate a bit for if triangle coincide with processor boundaries
            for (auto& bb : bbs[proci])
            {
                bb.grow(mergeDist_);
            }
            Pstream::listCombineGather(bbs[proci], plusEqOp<boundBox>());
        }
        Pstream::broadcast(bbs);
    }
    else if (distType_ == DISTRIBUTED)
    {
        // Fully distributed decomposition. Does not filter duplicate
        // triangles.
        if (!decomposer_().parallelAware())
        {
            FatalErrorInFunction
                << "The decomposition method " << decomposer_().typeName
                << " does not decompose in parallel."
                << " Please choose one that does." << exit(FatalError);
        }

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
                << " determining decomposition for " << s.size()
                << " centroids of surface " << searchableSurface::name()
                << endl;
        }

        // Triangle centres
        pointField triCentres(s.size());
        forAll(s, trii)
        {
            triCentres[trii] = s[trii].centre(s.points());
        }

        labelList distribution = decomposer_().decompose(triCentres);

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
                << " determining processor bounding boxes for "
                << searchableSurface::name() << endl;
        }

        // Find bounding box for all triangles on new distribution.
        forAll(s, trii)
        {
            treeBoundBox& bb = bbs[distribution[trii]][0];
            bb.add(s.points(), s[trii]);
        }

        // Now combine for all processors and convert to correct format.
        forAll(bbs, proci)
        {
            // Inflate a bit for if triangle coincide with processor boundaries
            for (auto& bb : bbs[proci])
            {
                bb.grow(mergeDist_);
            }
            Pstream::listCombineGather(bbs[proci], plusEqOp<boundBox>());
        }
        Pstream::broadcast(bbs);
    }
//    //// Tbd. initial duplicate filtering of border points only
//    if (distType_ == DISTRIBUTED)
//    {
//        // Fully distributed decomposition. Does not filter duplicate
//        // triangles.
//        if (!decomposer_().parallelAware())
//        {
//            FatalErrorInFunction
//                << "The decomposition method " << decomposer_().typeName
//                << " does not decompose in parallel."
//                << " Please choose one that does." << exit(FatalError);
//        }
//
//        if (debug)
//        {
//            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
//                << " determining decomposition for " << s.size()
//                << " centroids" << endl;
//        }
//        const pointField& points = s.points();
//
//        pointField triCentres(s.size());
//        forAll(s, trii)
//        {
//            triCentres[trii] = s[trii].centre(points);
//        }
//
//        // Collect all triangles not fully inside the current bb
//        DynamicList<label> borderTris(s.size()/Pstream::nProcs());
//
//        const List<treeBoundBox>& myBbs = procBb_[Pstream::myProcNo];
//
//        boolList includedFace;
//        overlappingSurface(s, myBbs, includedFace);
//        boolList internalOrBorderFace(includedFace);
//        forAll(s, trii)
//        {
//            if (includedFace[trii])
//            {
//                // Triangle is either inside or part-inside. Exclude fully
//                // inside triangles.
//                const labelledTri& f = s[trii];
//                const point& p0 = points[f[0]];
//                const point& p1 = points[f[1]];
//                const point& p2 = points[f[2]];
//                if
//                (
//                   !contains(myBbs, p0)
//                || !contains(myBbs, p1)
//                || !contains(myBbs, p2)
//                )
//                {
//                    borderTris.append(trii);
//                }
//            }
//        }
//
//        const label nBorderTris = borderTris.size();
//
//        Pout<< "Have " << borderTris.size() << " border triangles out of "
//            << s.size() << endl;
//
//        labelListList sendMap(Pstream::nProcs());
//        sendMap[0] = std::move(borderTris);
//
//        const mapDistribute map(std::move(sendMap));
//
//        // Gather all borderTris
//        //globalIndex globalBorderTris(borderTris.size());
//        //pointField globalBorderCentres(allCentres, borderTris);
//        //globalBorderTris.gather
//        //(
//        //    UPstream::worldComm,
//        //    UPstream::allProcs(UPstream::worldComm),
//        //    globalBorderCentres
//        //);
//        pointField globalBorderCentres(allCentres);
//        map.distribute(globalBorderCentres);
//
//        // Merge on master
//        labelList masterBorder(borderTris.size(), -1);
//        if (Pstream::master())
//        {
//            labelList allToMerged;
//            label nMerged = mergePoints
//            (
//                globalBorderCentres,
//                mergeDist_,
//                false,          // verbose = false
//                allToMerged
//            );
//
//            if (debug)
//            {
//                Pout<< "distributedTriSurfaceMesh::"
//                    << "independentlyDistributedBbs :"
//                    << " merged " << globalBorderCentres.size()
//                    << " border centroids down to " << nMerged << endl;
//            }
//
//            labelList mergedMaster(nMerged, -1);
//            isMaster.setSize(globalBorderCentres.size(), false);
//            forAll(allToMerged, i)
//            {
//                label mergedi = allToMerged[i];
//                if (mergedMaster[mergedi] == -1)
//                {
//                    mergedMaster[mergedi] = i;
//                    isMaster[i] = true;
//                }
//            }
//            forAll(allToMerged, i)
//            {
//                label mergedi = allToMerged[i];
//                masterBorder[i] = mergedMaster[mergedi];
//            }
//        }
//        //globalBorderTris.scatter
//        //(
//        //    UPstream::worldComm,
//        //    UPstream::allProcs(UPstream::worldComm),
//        //    isMasterPoint
//        //);
//        //boolList isMasterBorder(s.size(), false);
//        //forAll(borderTris, i)
//        //{
//        //    isMasterBorder[borderTris[i]] = isMasterPoint[i];
//        //}
//        map.reverseDistribute(s.size(), isMaster);
//
//        // Collect all non-border or master-border points
//        DynamicList<label> triMap(s.size());
//        forAll(s, trii)
//        {
//            if (includedFace[trii])
//            {
//                // Triangle is either inside or part-inside. Exclude fully
//                // inside triangles.
//                const labelledTri& f = s[trii];
//                const point& p0 = points[f[0]];
//                const point& p1 = points[f[1]];
//                const point& p2 = points[f[2]];
//
//                if
//                (
//                   contains(myBbs, p0)
//                && contains(myBbs, p1)
//                && contains(myBbs, p2)
//                )
//                {
//                    // Internal
//                    triMap.append(trii);
//                }
//                else if (isMasterBorder[trii])
//                {
//                    // Part overlapping and master of overlap
//                    triMap.append(trii);
//                }
//            }
//        }
//
//        pointField masterCentres(allCentres, triMap);
//
//        // Do the actual decomposition
//        labelList masterDistribution(decomposer_().decompose(masterCentres));
//
//        // Make map to get the decomposition from the master of each border
//        labelList borderGlobalMaster(borderTris.size());
//        forAll(borderTris, borderi)
//        {
//            borderGlobalMaster[borderi] = ..masterTri
//        }
//        mapDistribute map(globalBorderTris, borderGlobalMaster
//
//        // Send decomposition
//
//
//        if (debug)
//        {
//            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
//                << " determining processor bounding boxes" << endl;
//        }
//
//        // Find bounding box for all triangles on new distribution.
//        forAll(s, trii)
//        {
//            treeBoundBox& bb = bbs[distribution[trii]][0];
//            bb.add(s.points(), s[trii]);
//        }
//
//        // Now combine for all processors and convert to correct format.
//        forAll(bbs, proci)
//        {
//            Pstream::listCombineGather(bbs[proci], plusEqOp<boundBox>());
//        }
//        Pstream::broadcast(bbs);
//    }
    else
    {
        // Master-only decomposition. Filters duplicate triangles so repeatable.

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
                << " collecting all centroids for surface "
                << searchableSurface::name() << endl;
        }

        // Collect all triangle centres
        pointField allCentres(s.size());
        {
            forAll(s, trii)
            {
                allCentres[trii] = s[trii].centre(s.points());
            }
            globalTris().gather
            (
                UPstream::worldComm,
                UPstream::allProcs(UPstream::worldComm),
                allCentres
            );
        }

        // Determine local decomposition
        labelList distribution(s.size());
        {
            labelList allDistribution;
            if (Pstream::master())
            {
                labelList allToMerged;
                label nMerged = mergePoints
                (
                    allCentres,
                    mergeDist_,
                    false,          // verbose = false
                    allToMerged
                );

                if (debug)
                {
                    Pout<< "distributedTriSurfaceMesh::"
                        << "independentlyDistributedBbs :"
                        << " surface:" << searchableSurface::name()
                        << " merged " << allCentres.size()
                        << " centroids down to " << nMerged << endl;
                }

                // Collect merged points
                pointField mergedPoints(nMerged);
                UIndirectList<point>(mergedPoints, allToMerged) = allCentres;

                // Decompose merged centres
                const bool oldParRun = UPstream::parRun(false);
                labelList mergedDist(decomposer_().decompose(mergedPoints));
                UPstream::parRun(oldParRun);  // Restore parallel state

                // Convert to all
                allDistribution = UIndirectList<label>
                (
                    mergedDist,
                    allToMerged
                );
            }

            // Scatter back to processors
            globalTris().scatter
            (
                UPstream::worldComm,
                UPstream::allProcs(UPstream::worldComm),
                allDistribution,
                distribution
            );
            if (debug)
            {
                Pout<< "distributedTriSurfaceMesh::"
                    << "independentlyDistributedBbs :"
                    << " determined decomposition" << endl;
            }
        }

        // Find bounding box for all triangles on new distribution.
        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::independentlyDistributedBbs :"
                << " determining processor bounding boxes for "
                << searchableSurface::name() << endl;
        }

        forAll(s, trii)
        {
            treeBoundBox& bb = bbs[distribution[trii]][0];
            bb.add(s.points(), s[trii]);
        }

        // Now combine for all processors and convert to correct format.
        forAll(bbs, proci)
        {
            // Inflate a bit for if triangle coincide with processor boundaries
            for (auto& bb : bbs[proci])
            {
                bb.grow(mergeDist_);
            }
            Pstream::listCombineGather(bbs[proci], plusEqOp<boundBox>());
        }
        Pstream::broadcast(bbs);
    }
    return bbs;
}


// Does any part of triangle overlap bb.
bool Foam::distributedTriSurfaceMesh::overlaps
(
    const List<treeBoundBox>& bbs,
    const triPointRef& tri
)
{
    treeBoundBox triBb(tri.a());
    triBb.add(tri.b());
    triBb.add(tri.c());

    for (const treeBoundBox& bb : bbs)
    {
        // Exact test of triangle intersecting bb

        // Quick rejection. If whole bounding box of tri is outside cubeBb then
        // there will be no intersection.

        if (bb.overlaps(triBb) && bb.intersects(tri))
        {
            return true;
        }
    }
    return false;
}


void Foam::distributedTriSurfaceMesh::subsetMeshMap
(
    const triSurface& s,
    const boolList& include,
    const label nIncluded,
    labelList& newToOldPoints,
    labelList& oldToNewPoints,
    labelList& newToOldFaces
)
{
    newToOldFaces.setSize(nIncluded);
    newToOldPoints.setSize(s.points().size());
    oldToNewPoints.setSize(s.points().size());
    oldToNewPoints = -1;
    {
        label facei = 0;
        label pointi = 0;

        forAll(include, oldFacei)
        {
            if (include[oldFacei])
            {
                // Store new faces compact
                newToOldFaces[facei++] = oldFacei;

                // Renumber labels for face
                for (const label oldPointi : s[oldFacei])
                {
                    if (oldToNewPoints[oldPointi] == -1)
                    {
                        oldToNewPoints[oldPointi] = pointi;
                        newToOldPoints[pointi++] = oldPointi;
                    }
                }
            }
        }
        newToOldPoints.setSize(pointi);
    }
}


Foam::triSurface Foam::distributedTriSurfaceMesh::subsetMesh
(
    const triSurface& s,
    const labelList& newToOldPoints,
    const labelList& oldToNewPoints,
    const labelList& newToOldFaces
)
{
    // Extract points
    pointField newPoints(newToOldPoints.size());
    forAll(newToOldPoints, i)
    {
        newPoints[i] = s.points()[newToOldPoints[i]];
    }
    // Extract faces
    List<labelledTri> newTriangles(newToOldFaces.size());

    forAll(newToOldFaces, i)
    {
        // Get old vertex labels
        const labelledTri& tri = s[newToOldFaces[i]];

        newTriangles[i][0] = oldToNewPoints[tri[0]];
        newTriangles[i][1] = oldToNewPoints[tri[1]];
        newTriangles[i][2] = oldToNewPoints[tri[2]];
        newTriangles[i].region() = tri.region();
    }

    // Reuse storage.
    return triSurface(newTriangles, s.patches(), newPoints, true);
}


Foam::triSurface Foam::distributedTriSurfaceMesh::subsetMesh
(
    const triSurface& s,
    const boolList& include,
    labelList& newToOldPoints,
    labelList& newToOldFaces
)
{
    label n = 0;

    forAll(include, i)
    {
        if (include[i])
        {
            n++;
        }
    }

    labelList oldToNewPoints;
    subsetMeshMap
    (
        s,
        include,
        n,
        newToOldPoints,
        oldToNewPoints,
        newToOldFaces
    );

    return subsetMesh
    (
        s,
        newToOldPoints,
        oldToNewPoints,
        newToOldFaces
    );
}


Foam::triSurface Foam::distributedTriSurfaceMesh::subsetMesh
(
    const triSurface& s,
    const labelList& newToOldFaces,
    labelList& newToOldPoints
)
{
    const boolList include
    (
        ListOps::createWithValue<bool>(s.size(), newToOldFaces, true, false)
    );

    newToOldPoints.setSize(s.points().size());
    labelList oldToNewPoints(s.points().size(), -1);
    {
        label pointi = 0;

        forAll(include, oldFacei)
        {
            if (include[oldFacei])
            {
                // Renumber labels for face
                for (const label oldPointi : s[oldFacei])
                {
                    if (oldToNewPoints[oldPointi] == -1)
                    {
                        oldToNewPoints[oldPointi] = pointi;
                        newToOldPoints[pointi++] = oldPointi;
                    }
                }
            }
        }
        newToOldPoints.setSize(pointi);
    }

    return subsetMesh
    (
        s,
        newToOldPoints,
        oldToNewPoints,
        newToOldFaces
    );
}


Foam::label Foam::distributedTriSurfaceMesh::findTriangle
(
    const List<labelledTri>& allFaces,
    const labelListList& allPointFaces,
    const labelledTri& otherF
)
{
    // allFaces connected to otherF[0]
    const labelList& pFaces = allPointFaces[otherF[0]];

    forAll(pFaces, i)
    {
        const labelledTri& f = allFaces[pFaces[i]];

        if (f.region() == otherF.region())
        {
            // Find index of otherF[0]
            label fp0 = f.find(otherF[0]);
            // Check rest of triangle in same order
            label fp1 = f.fcIndex(fp0);
            label fp2 = f.fcIndex(fp1);

            if (f[fp1] == otherF[1] && f[fp2] == otherF[2])
            {
                return pFaces[i];
            }
        }
    }
    return -1;
}


// Merge into allSurf.
void Foam::distributedTriSurfaceMesh::merge
(
    const scalar mergeDist,
    const List<labelledTri>& subTris,
    const pointField& subPoints,

    List<labelledTri>& allTris,
    pointField& allPoints,

    labelList& faceConstructMap,
    labelList& pointConstructMap
)
{
    labelList subToAll;
    matchPoints
    (
        subPoints,
        allPoints,
        scalarField(subPoints.size(), mergeDist),   // match distance
        false,                                      // verbose
        pointConstructMap
    );

    label nOldAllPoints = allPoints.size();


    // Add all unmatched points
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    label allPointi = nOldAllPoints;
    forAll(pointConstructMap, pointi)
    {
        if (pointConstructMap[pointi] == -1)
        {
            pointConstructMap[pointi] = allPointi++;
        }
    }

    if (allPointi > nOldAllPoints)
    {
        allPoints.setSize(allPointi);

        forAll(pointConstructMap, pointi)
        {
            if (pointConstructMap[pointi] >= nOldAllPoints)
            {
                allPoints[pointConstructMap[pointi]] = subPoints[pointi];
            }
        }
    }


    // To check whether triangles are same we use pointFaces.
    labelListList allPointFaces;
    invertManyToMany(nOldAllPoints, allTris, allPointFaces);


    // Add all unmatched triangles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label allTrii = allTris.size();
    allTris.setSize(allTrii + subTris.size());

    faceConstructMap.setSize(subTris.size());

    forAll(subTris, trii)
    {
        const labelledTri& subTri = subTris[trii];

        // Get triangle in new numbering
        labelledTri mappedTri
        (
            pointConstructMap[subTri[0]],
            pointConstructMap[subTri[1]],
            pointConstructMap[subTri[2]],
            subTri.region()
        );


        // Check if all points were already in surface
        bool fullMatch = true;

        forAll(mappedTri, fp)
        {
            if (mappedTri[fp] >= nOldAllPoints)
            {
                fullMatch = false;
                break;
            }
        }

        if (fullMatch)
        {
            // All three points are mapped to old points. See if same
            // triangle.
            label i = findTriangle
            (
                allTris,
                allPointFaces,
                mappedTri
            );

            if (i == -1)
            {
                // Add
                faceConstructMap[trii] = allTrii;
                allTris[allTrii] = mappedTri;
                allTrii++;
            }
            else
            {
                faceConstructMap[trii] = i;
            }
        }
        else
        {
            // Add
            faceConstructMap[trii] = allTrii;
            allTris[allTrii] = mappedTri;
            allTrii++;
        }
    }
    allTris.setSize(allTrii);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh
(
    const IOobject& io,
    const triSurface& s,
    const dictionary& dict
)
:
    triSurfaceMesh(io, s),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            searchableSurface::NO_READ,
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        ),
        dict
    ),
    currentDistType_(FROZEN)    // only used to trigger re-distribution
{
    // Read from the provided dictionary
    read();

    bounds().reduce();

    if (debug)
    {
        InfoInFunction << "Constructed from triSurface "
            << searchableSurface::name()
            << endl;
        writeStats(Info);

        labelList nTris
        (
            UPstream::listGatherValues<label>(triSurface::size())
        );

        if (Pstream::master())
        {
            Info<< endl<< "\tproc\ttris\tbb" << endl;
            forAll(nTris, proci)
            {
                Info<< '\t' << proci << '\t' << nTris[proci]
                    << '\t' << procBb_[proci] << endl;
            }
            Info<< endl;
        }
    }
}


Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh(const IOobject& io)
:
    triSurfaceMesh
    (
        IOobject
        (
            io.name(),
            findLocalInstance(io),  // findInstance with parent searching
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject()
        ),
        triSurfaceMesh::masterOnly    // allow parent searching
    ),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            IOobjectOption::lazierRead(searchableSurface::readOpt()),
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        ),
        dictionary()
    ),
    currentDistType_(FROZEN)    // only used to trigger re-distribution
{
    // Read from the local, decomposed dictionary
    read();

    bounds().reduce();

    // Try and find out where the triSurface was loaded from. On slave this
    // might give empty filename.
    const fileName actualFile(triSurfaceMesh::findFile(io, true));

    if
    (
        (
            actualFile.empty()
         || actualFile != io.localFilePath(triSurfaceMesh::typeName)
        )
     && (distType_ == INDEPENDENT || distType_ == DISTRIBUTED)
    )
    {
        DebugInFunction
            << "Read distributedTriSurface " << io.name()
            << " from parent path " << actualFile << endl;

        if (Pstream::parRun())
        {
            // Distribute (checks that distType != currentDistType_ so should
            // always trigger re-distribution)
            List<treeBoundBox> bbs;
            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            distribute
            (
                bbs,
                true,       // keep unmapped triangles
                faceMap,
                pointMap
            );
        }
    }
    else
    {
        if (debug)
        {
            InfoInFunction
                << "Read distributedTriSurface " << io.name()
                << " from actual path " << actualFile << ':' << endl;

            labelList nTris
            (
                UPstream::listGatherValues<label>(triSurface::size())
            );

            if (Pstream::master())
            {
                Info<< endl<< "\tproc\ttris\tbb" << endl;
                forAll(nTris, proci)
                {
                    Info<< '\t' << proci << '\t' << nTris[proci]
                        << '\t' << procBb_[proci] << endl;
                }
                Info<< endl;
            }
        }
    }
    if (debug)
    {
        InfoInFunction
            << "Read distributedTriSurface " << io.name() << ':' << endl;
        writeStats(Info);
    }
}


Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh
    (
        IOobject
        (
            io.name(),
            findLocalInstance(io),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject()
        ),
        dict,
        triSurfaceMesh::masterOnly
    ),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            IOobjectOption::lazierRead(searchableSurface::readOpt()),
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        ),
        dictionary()
    ),
    currentDistType_(FROZEN)    // only used to trigger re-distribution
{
    // Read from the local, decomposed dictionary
    read();

    // Optionally override settings from provided dictionary
    {
        // Wanted distribution type
        distributionTypeNames_.readIfPresent
        (
            "distributionType",
            dict,
            distType_
        );

        // Merge distance
        dict.readIfPresent("mergeDistance", mergeDist_);

        // Distribution type
        bool closed;
        if (dict.readIfPresent<bool>("closed", closed))
        {
            surfaceClosed_ = closed;
        }

        outsideVolType_ = volumeType::names.getOrDefault
        (
            "outsideVolumeType",
            dict,
            outsideVolType_
        );
    }


    bounds().reduce();

    // Try and find out where the triSurface was loaded from. On slave this
    // might give empty filename.
    const fileName actualFile(triSurfaceMesh::findFile(io, dict, true));

    if
    (
        (
            actualFile.empty()
         || actualFile != io.localFilePath(triSurfaceMesh::typeName)
        )
     && (distType_ == INDEPENDENT || distType_ == DISTRIBUTED)
    )
    {
        DebugInFunction
            << "Read distributedTriSurface " << io.name()
            << " from parent path " << actualFile
            << " and dictionary" << endl;

        if (Pstream::parRun())
        {
            // Distribute (checks that distType != currentDistType_ so should
            // always trigger re-distribution)
            List<treeBoundBox> bbs;
            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            distribute
            (
                bbs,
                true,       // keep unmapped triangles
                faceMap,
                pointMap
            );
        }
    }
    else
    {
        if (debug)
        {
            InfoInFunction
                << "Read distributedTriSurface " << io.name()
                << " from actual path " << actualFile
                << " and dictionary:" << endl;

            labelList nTris
            (
                UPstream::listGatherValues<label>(triSurface::size())
            );

            if (Pstream::master())
            {
                Info<< endl<< "\tproc\ttris\tbb" << endl;
                forAll(nTris, proci)
                {
                    Info<< '\t' << proci << '\t' << nTris[proci]
                        << '\t' << procBb_[proci] << endl;
                }
                Info<< endl;
            }
        }
    }
    if (debug)
    {
        InfoInFunction
            << "Read distributedTriSurface " << io.name() << ':' << endl;
        writeStats(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributedTriSurfaceMesh::~distributedTriSurfaceMesh()
{
    clearOut();
}


void Foam::distributedTriSurfaceMesh::clearOut()
{
    globalTris_.clear();
    triSurfaceMesh::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalIndex& Foam::distributedTriSurfaceMesh::globalTris() const
{
    if (!globalTris_)
    {
        globalTris_.reset(new globalIndex(triSurface::size()));
    }
    return *globalTris_;
}


//void Foam::distributedTriSurfaceMesh::findNearest
//(
//    const pointField& samples,
//    const scalarField& nearestDistSqr,
//    List<pointIndexHit>& info
//) const
//{
//    if (!Pstream::parRun())
//    {
//        triSurfaceMesh::findNearest(samples, nearestDistSqr, info);
//        return;
//    }
//
//    addProfiling
//    (
//        findNearest,
//        "distributedTriSurfaceMesh::findNearest"
//    );
//
//    if (debug)
//    {
//        Pout<< "distributedTriSurfaceMesh::findNearest :"
//            << " trying to find nearest for " << samples.size()
//            << " samples with max sphere "
//            << (samples.size() ? Foam::sqrt(max(nearestDistSqr)) : Zero)
//            << endl;
//    }
//
//
//    const indexedOctree<treeDataTriSurface>& octree = tree();
//
//    // Important:force synchronised construction of indexing
//    const globalIndex& triIndexer = globalTris();
//
//
//    // Initialise
//    // ~~~~~~~~~~
//
//    info.setSize(samples.size());
//    forAll(info, i)
//    {
//        info[i].setMiss();
//    }
//
//
//
//    // Do any local queries
//    // ~~~~~~~~~~~~~~~~~~~~
//
//    label nLocal = 0;
//
//    {
//        // Work array - whether processor bb overlaps the bounding sphere.
//        boolList procBbOverlaps(Pstream::nProcs());
//
//        forAll(samples, i)
//        {
//            // Find the processor this sample+radius overlaps.
//            label nProcs = calcOverlappingProcs
//            (
//                samples[i],
//                nearestDistSqr[i],
//                procBbOverlaps
//            );
//
//            // Overlaps local processor?
//            if (procBbOverlaps[Pstream::myProcNo()])
//            {
//                info[i] = octree.findNearest(samples[i], nearestDistSqr[i]);
//                if (info[i].hit())
//                {
//                    if
//                    (
//                        surfaceClosed_
//                    && !contains(procBb_[proci], info[i].point())
//                    )
//                    {
//                        // Nearest point is not on local processor so the
//                        // the triangle is only there because some other bit
//                        // of it
//                        // is on it. Assume there is another processor that
//                        // holds
//                        // the full surrounding of the triangle so we can
//                        // clear  this particular nearest.
//                        info[i].setMiss();
//                        info[i].setIndex(-1);
//                    }
//                    else
//                    {
//                        info[i].setIndex
//                        (triIndexer.toGlobal(info[i].index()));
//                    }
//                }
//                if (nProcs == 1)
//                {
//                    // Fully local
//                    nLocal++;
//                }
//            }
//        }
//    }
//
//
//    if
//    (
//        Pstream::parRun()
//     && (
//            returnReduce(nLocal, sumOp<label>())
//          < returnReduce(samples.size(), sumOp<label>())
//        )
//    )
//    {
//        // Not all can be resolved locally. Build queries and map, send over
//        // queries, do intersections, send back and merge.
//
//        // Calculate queries and exchange map
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//        pointField allCentres;
//        scalarField allRadiusSqr;
//        labelList allSegmentMap;
//        autoPtr<mapDistribute> mapPtr
//        (
//            calcLocalQueries
//            (
//                false,    // exclude local processor since already done above
//                samples,
//                nearestDistSqr,
//
//                allCentres,
//                allRadiusSqr,
//                allSegmentMap
//            )
//        );
//        const mapDistribute& map = mapPtr();
//
//
//        // swap samples to local processor
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//        map.distribute(allCentres);
//        map.distribute(allRadiusSqr);
//
//
//        // Do my tests
//        // ~~~~~~~~~~~
//
//        List<pointIndexHit> allInfo(allCentres.size());
//        forAll(allInfo, i)
//        {
//            allInfo[i] = octree.findNearest
//            (
//                allCentres[i],
//                allRadiusSqr[i]
//            );
//            if (allInfo[i].hit())
//            {
//                // We don't know if the nearest is on an edge/point. If
//                // this is the case we preferentially want to return the
//                // index on the processor that holds all surrounding triangles
//                // so we can do e.g. follow-on inside/outside tests
//                if
//                (
//                    surfaceClosed_
//                && !contains(procBb_[Pstream::myProcNo()], allInfo[i].point())
//                )
//                {
//                    // Nearest point is not on local processor so the
//                    // the triangle is only there because some other bit of it
//                    // is on it. Assume there is another processor that holds
//                    // the full surrounding of the triangle so we can clear
//                    // this particular nearest.
//                    allInfo[i].setMiss();
//                    allInfo[i].setIndex(-1);
//                }
//                else
//                {
//                    allInfo[i].setIndex
//                    (
//                        triIndexer.toGlobal(allInfo[i].index())
//                    );
//                }
//            }
//        }
//
//
//        // Send back results
//        // ~~~~~~~~~~~~~~~~~
//
//        map.reverseDistribute(allSegmentMap.size(), allInfo);
//
//
//        // Extract information
//        // ~~~~~~~~~~~~~~~~~~~
//
//        forAll(allInfo, i)
//        {
//            if (allInfo[i].hit())
//            {
//                label pointi = allSegmentMap[i];
//
//                if (!info[pointi].hit())
//                {
//                    // No intersection yet so take this one
//                    info[pointi] = allInfo[i];
//                }
//                else
//                {
//                    // Nearest intersection
//                    if
//                    (
//                        samples[pointi].distSqr(allInfo[i].point())
//                      < samples[pointi].distSqr(info[pointi].point())
//                    )
//                    {
//                        info[pointi] = allInfo[i];
//                    }
//                }
//            }
//        }
//    }
//}


void Foam::distributedTriSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::findNearest(samples, nearestDistSqr, info);
        return;
    }

    addProfiling
    (
        findNearest,
        "distributedTriSurfaceMesh::findNearest"
    );

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findNearest :"
            << " surface " << searchableSurface::name()
            << " trying to find nearest for " << samples.size()
            << " samples with max sphere "
            << (samples.size() ? Foam::sqrt(max(nearestDistSqr)) : Zero)
            << endl;
    }

    const globalIndex& triIndexer = globalTris();

    // Two-pass searching:
    // 1. send the sample to the processor whose bb contains it. This is
    //    most likely also the one that holds the nearest triangle. (In case
    //    there is no containing processor send to nearest processors. Note
    //    that this might cause a lot of traffic if this is likely)
    //    Send the resulting nearest point back.
    // 2. with the find from 1 look at which other processors might have a
    //    better triangle. Since hopefully step 1) will have produced a tight
    //    bounding box this should limit the amount of points to be retested


    // 1. Test samples on processor(s) that contains them
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<mapDistribute> map1Ptr;
    scalarField distSqr(nearestDistSqr);
    boolList procContains(Pstream::nProcs(), false);
    boolList procOverlaps(Pstream::nProcs(), false);

    label nOutside = 0;
    {
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());
        // Pre-size. Assume samples are uniformly distributed
        forAll(dynSendMap, proci)
        {
            dynSendMap[proci].reserve(samples.size()/Pstream::nProcs());
        }

        forAll(samples, samplei)
        {
            label minProci = -1;
            Tuple2<label, scalar> best = findBestProcs
            (
                samples[samplei],
                distSqr[samplei],
                procContains,
                procOverlaps,
                minProci
            );

            label nContains = 0;
            forAll(procBb_, proci)
            {
                if (procContains[proci])
                {
                    nContains++;
                    dynSendMap[proci].append(samplei);
                    distSqr[samplei] = best.second();
                }
            }
            if (nContains == 0)
            {
                nOutside++;
                // Sample is outside all bb. Choices:
                //  - send to all processors
                //  - send to single processor

                //forAll(procOverlaps[proci])
                //{
                //    if (procOverlaps[proci])
                //    {
                //        dynSendMap[proci].append(samplei);
                //        distSqr[samplei] = best.second();
                //    }
                //}
                if (minProci != -1)
                {
                    dynSendMap[minProci].append(samplei);
                    distSqr[samplei] = best.second();
                }
            }
        }

        labelListList sendMap(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }
        map1Ptr.reset(new mapDistribute(std::move(sendMap)));
    }
    const mapDistribute& map1 = *map1Ptr;


    if (debug)
    {
        Pout<< searchableSurface::name() << " Pass1:"
            << " of " << samples.size() << " samples sending to" << endl;
        label nSend = 0;
        forAll(map1.subMap(), proci)
        {
            Pout<< "    " << proci << "\t" << map1.subMap()[proci].size()
                << endl;
            nSend += map1.subMap()[proci].size();
        }
        Pout<< "    sum\t" << nSend << endl
            << "    outside\t" << nOutside << endl;
    }


    List<nearestAndDist> nearestInfo;
    {
        // Get the points I need to test and test locally
        pointField localPoints(samples);
        map1.distribute(localPoints);
        scalarField localDistSqr(distSqr);
        map1.distribute(localDistSqr);
        List<pointIndexHit> localInfo;
        triSurfaceMesh::findNearest(localPoints, localDistSqr, localInfo);
        convertTriIndices(localInfo);

        // Pack into structure for combining information from multiple
        // processors
        nearestInfo.setSize(localInfo.size());
        nearestInfo = nearestAndDist(pointIndexHit(), Foam::sqr(GREAT));

        label nHit = 0;
        label nIgnoredHit = 0;

        forAll(nearestInfo, i)
        {
            const pointIndexHit& info = localInfo[i];
            if (info.hit())
            {
                nHit++;

                if
                (
                    surfaceClosed_
                && !contains(procBb_[Pstream::myProcNo()], info.point())
                )
                {
                    // Nearest point is not on local processor so the
                    // the triangle is only there because some other bit
                    // of it is on it. Assume there is another processor that
                    // holds the full surrounding of the triangle so we can
                    // ignore this particular nearest.
                    nIgnoredHit++;
                }
                else
                {
                    nearestAndDist& ni = nearestInfo[i];
                    ni.first() = info;
                    ni.second() = info.point().distSqr(localPoints[i]);
                }
            }
        }

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::findNearest :"
                << " surface " << searchableSurface::name()
                << " searched locally for " << localPoints.size()
                << " samples with max sphere "
                << (localDistSqr.size() ? Foam::sqrt(max(localDistSqr)) : Zero)
                << " found hits:" << nHit
                << " of which outside local bb:" << nIgnoredHit
                << endl;
        }
    }

    // Send back to originating processor. Choose best if sent to multiple
    // processors. Note that afterwards all unused entries have the unique
    // value nearestZero (distance < 0). This is used later on to see if
    // the sample was sent to any processor.
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        samples.size(),
        map1.constructMap(),
        map1.constructHasFlip(),
        map1.subMap(),
        map1.subHasFlip(),
        nearestInfo,
        nearestZero,
        nearestEqOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map1.comm()
    );


    // 2. Test samples on other processor(s) that overlap
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Now we have (in nearestInfo) for every input sample the current best
    // hit (on the processor that originates the sample). See if we can
    // improve it by sending the queries to any other processors

    autoPtr<mapDistribute> map2Ptr;
    {
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());

        // Work array - whether processor bb overlaps the bounding sphere.
        boolList procBbOverlaps(Pstream::nProcs());

        // label nFound = 0;

        forAll(nearestInfo, samplei)
        {
            const point& sample = samples[samplei];
            const nearestAndDist& ni = nearestInfo[samplei];
            const pointIndexHit& info = ni.first();

            // if (info.hit())
            // {
            //     nFound++;
            // }

            scalar d2 =
            (
                info.hit()
              ? ni.second()
              : distSqr[samplei]
            );

            label hitProci =
            (
                info.hit()
              ? triIndexer.whichProcID(info.index())
              : -1
            );

            // Find the processors this sample+radius overlaps.
            calcOverlappingProcs(sample, d2, procBbOverlaps);

            forAll(procBbOverlaps, proci)
            {
                if (procBbOverlaps[proci])
                {
                    // Check this sample wasn't already handled above. This
                    // could be improved since the sample might have been
                    // searched on multiple processors. We now only exclude the
                    // processor where the point was inside.
                    if (proci != hitProci)
                    {
                        dynSendMap[proci].append(samplei);
                    }
                }
            }
        }

        labelListList sendMap(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }
        map2Ptr.reset(new mapDistribute(std::move(sendMap)));
    }

    const mapDistribute& map2 = map2Ptr();

    if (debug)
    {
        Pout << searchableSurface::name() << " Pass2:"
            << " of " << samples.size() << " samples sending to" << endl;
        label nSend = 0;
        forAll(map2.subMap(), proci)
        {
            Pout<< "    " << proci << "\t" << map2.subMap()[proci].size()
                << endl;
            nSend += map2.subMap()[proci].size();
        }
        Pout<< "    sum\t" << nSend << endl;
    }

    // Send samples and current best distance
    pointField localSamples(samples);
    map2.distribute(localSamples);
    scalarField localDistSqr(distSqr);
    forAll(nearestInfo, samplei)
    {
        const nearestAndDist& ni = nearestInfo[samplei];
        if (ni.first().hit())
        {
            localDistSqr[samplei] = ni.second();
        }
    }
    map2.distribute(localDistSqr);

    // Do local test
    List<pointIndexHit> localInfo;
    triSurfaceMesh::findNearest(localSamples, localDistSqr, localInfo);
    convertTriIndices(localInfo);

    // Pack and send back
    List<nearestAndDist> localBest(localSamples.size());
    label nHit = 0;
    label nIgnoredHit = 0;
    forAll(localInfo, i)
    {
        const pointIndexHit& info = localInfo[i];
        if (info.hit())
        {
            nHit++;
            if
            (
                surfaceClosed_
            && !contains(procBb_[Pstream::myProcNo()], info.point())
            )
            {
                // See above
                nIgnoredHit++;
            }
            else
            {
                nearestAndDist& ni = localBest[i];
                ni.first() = info;
                ni.second() = info.point().distSqr(localSamples[i]);
            }
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findNearest :"
            << " surface " << searchableSurface::name()
            << " searched locally for " << localSamples.size()
            << " samples with max sphere "
            << (localDistSqr.size() ? Foam::sqrt(max(localDistSqr)) : Zero)
            << " found hits:" << nHit
            << " of which outside local bb:" << nIgnoredHit
            << endl;
    }

    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        samples.size(),
        map2.constructMap(),
        map2.constructHasFlip(),
        map2.subMap(),
        map2.subHasFlip(),
        localBest,
        nearestZero,
        nearestEqOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map2.comm()
    );

    // Combine with nearestInfo
    info.setSize(samples.size());
    forAll(samples, samplei)
    {
        nearestAndDist ni(nearestInfo[samplei]);
        nearestEqOp()(ni, localBest[samplei]);

        info[samplei] = ni.first();
    }
}


void Foam::distributedTriSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const labelList& regionIndices,
    List<pointIndexHit>& info
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::findNearest
        (
            samples,
            nearestDistSqr,
            regionIndices,
            info
        );
        return;
    }

    addProfiling
    (
        findNearestRegion,
        "distributedTriSurfaceMesh::findNearestRegion"
    );

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findNearest :"
            << " surface " << searchableSurface::name()
            << " trying to find nearest and region for " << samples.size()
            << " samples with max sphere "
            << (samples.size() ? Foam::sqrt(max(nearestDistSqr)) : Zero)
            << endl;
    }

    if (regionIndices.empty())
    {
        findNearest(samples, nearestDistSqr, info);
    }
    else
    {
        // Calculate queries and exchange map
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        pointField allCentres;
        scalarField allRadiusSqr;
        labelList allSegmentMap;
        autoPtr<mapDistribute> mapPtr
        (
            calcLocalQueries
            (
                true,      // also send to local processor
                samples,
                nearestDistSqr,

                allCentres,
                allRadiusSqr,
                allSegmentMap
            )
        );
        const mapDistribute& map = mapPtr();


        // swap samples to local processor
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        map.distribute(allCentres);
        map.distribute(allRadiusSqr);


        // Do my tests
        // ~~~~~~~~~~~

        List<pointIndexHit> allInfo(allCentres.size());
        triSurfaceMesh::findNearest
        (
            allCentres,
            allRadiusSqr,
            regionIndices,
            allInfo
        );
        convertTriIndices(allInfo);

        forAll(allInfo, i)
        {
            if (allInfo[i].hit())
            {
                if
                (
                    surfaceClosed_
                && !contains(procBb_[Pstream::myProcNo()], allInfo[i].point())
                )
                {
                    // Nearest point is not on local processor so the
                    // the triangle is only there because some other bit of it
                    // is on it. Assume there is another processor that holds
                    // the full surrounding of the triangle so we can clear
                    // this particular nearest.
                    allInfo[i].setMiss();
                    allInfo[i].setIndex(-1);
                }
            }
        }


        // Send back results
        // ~~~~~~~~~~~~~~~~~

        map.reverseDistribute(allSegmentMap.size(), allInfo);


        // Extract information
        // ~~~~~~~~~~~~~~~~~~~

        forAll(allInfo, i)
        {
            if (allInfo[i].hit())
            {
                label pointi = allSegmentMap[i];

                if (!info[pointi].hit())
                {
                    // No intersection yet so take this one
                    info[pointi] = allInfo[i];
                }
                else
                {
                    // Nearest intersection
                    if
                    (
                        samples[pointi].distSqr(allInfo[i].point())
                      < samples[pointi].distSqr(info[pointi].point())
                    )
                    {
                        info[pointi] = allInfo[i];
                    }
                }
            }
        }
    }
}


void Foam::distributedTriSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::findLine(start, end, info);
    }
    else
    {
        findLine
        (
            true,   // nearestIntersection
            start,
            end,
            info
        );
    }
}


void Foam::distributedTriSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::findLineAny(start, end, info);
    }
    else
    {
        findLine
        (
            true,   // nearestIntersection
            start,
            end,
            info
        );
    }
}


void Foam::distributedTriSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::findLineAll(start, end, info);
        return;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findLineAll :"
            << " surface " << searchableSurface::name()
            << " intersecting with "
            << start.size() << " rays" << endl;
    }

    addProfiling
    (
        findLineAll,
        "distributedTriSurfaceMesh::findLineAll"
    );

    // Reuse fineLine. We could modify all of findLine to do multiple
    // intersections but this would mean a lot of data transferred so
    // for now we just find nearest intersection and retest from that
    // intersection onwards.

    // Work array.
    List<pointIndexHit> hitInfo(start.size());

    findLine
    (
        true,   // nearestIntersection
        start,
        end,
        hitInfo
    );

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (SMALL is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress
    const vectorField dirVec(end-start);
    const scalarField magSqrDirVec(magSqr(dirVec));
    const vectorField smallVec
    (
        ROOTSMALL*dirVec
      + vector::uniform(ROOTVSMALL)
    );

    // Copy to input and compact any hits
    labelList pointMap(start.size());
    pointField e0(start.size());
    pointField e1(start.size());
    label compacti = 0;

    info.setSize(hitInfo.size());
    forAll(hitInfo, pointi)
    {
        if (hitInfo[pointi].hit())
        {
            info[pointi].setSize(1);
            info[pointi][0] = hitInfo[pointi];

            point pt = hitInfo[pointi].point() + smallVec[pointi];

            if (((pt-start[pointi])&dirVec[pointi]) <= magSqrDirVec[pointi])
            {
                e0[compacti] = pt;
                e1[compacti] = end[pointi];
                pointMap[compacti] = pointi;
                compacti++;
            }
        }
        else
        {
            info[pointi].clear();
        }
    }

    e0.setSize(compacti);
    e1.setSize(compacti);
    pointMap.setSize(compacti);


    label iter = 0;
    while (returnReduceOr(e0.size()))
    {
        findLine
        (
            true,   // nearestIntersection
            e0,
            e1,
            hitInfo
        );


        // Extract
        label compacti = 0;
        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                label pointi = pointMap[i];

                label sz = info[pointi].size();
                info[pointi].setSize(sz+1);
                info[pointi][sz] = hitInfo[i];

                point pt = hitInfo[i].point() + smallVec[pointi];

                // Check current coordinate along ray
                scalar d = ((pt-start[pointi])&dirVec[pointi]);

                // Note check for d>0. Very occasionally the octree will find
                // an intersection to the left of the ray due to tolerances.
                if (d > 0 && d <= magSqrDirVec[pointi])
                {
                    e0[compacti] = pt;
                    e1[compacti] = end[pointi];
                    pointMap[compacti] = pointi;
                    compacti++;
                }
            }
        }

        // Trim
        e0.setSize(compacti);
        e1.setSize(compacti);
        pointMap.setSize(compacti);

        iter++;

        if (iter == 1000)
        {
            Pout<< "distributedTriSurfaceMesh::findLineAll :"
                << " Exiting loop due to excessive number of"
                << " intersections along ray"
                << " start:" << UIndirectList<point>(start, pointMap)
                << " end:" << UIndirectList<point>(end, pointMap)
                << " e0:" << UIndirectList<point>(e0, pointMap)
                << " e1:" << UIndirectList<point>(e1, pointMap)
                << endl;
            break;
        }
    }
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findLineAll :"
            << " surface " << searchableSurface::name()
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::distributedTriSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getRegion :"
            << " surface " << searchableSurface::name()
            << " getting region for "
            << info.size() << " triangles" << endl;
    }

    addProfiling(getRegion, "distributedTriSurfaceMesh::getRegion");

    if (!Pstream::parRun())
    {
        region.setSize(info.size());
        forAll(info, i)
        {
            if (info[i].hit())
            {
                region[i] = triSurface::operator[](info[i].index()).region();
            }
            else
            {
                region[i] = -1;
            }
        }

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::getRegion :"
                << " surface " << searchableSurface::name()
                << " finished getting region for "
                << info.size() << " triangles" << endl;
        }

        return;
    }

    // Get query data (= local index of triangle)
    // ~~~~~~~~~~~~~~

    labelList triangleIndex(info.size());
    autoPtr<mapDistribute> mapPtr
    (
        localQueries
        (
            info,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();


    // Do my tests
    // ~~~~~~~~~~~

    const triSurface& s = static_cast<const triSurface&>(*this);

    region.setSize(triangleIndex.size());

    forAll(triangleIndex, i)
    {
        label trii = triangleIndex[i];
        region[i] = s[trii].region();
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    map.reverseDistribute(info.size(), region);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getRegion :"
            << " surface " << searchableSurface::name()
            << " finished getting region for "
            << info.size() << " triangles" << endl;
    }
}


void Foam::distributedTriSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::getNormal(info, normal);
        return;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getNormal :"
            << " surface " << searchableSurface::name()
            << " getting normal for "
            << info.size() << " triangles" << endl;
    }

    addProfiling(getNormal, "distributedTriSurfaceMesh::getNormal");


    // Get query data (= local index of triangle)
    // ~~~~~~~~~~~~~~

    labelList triangleIndex(info.size());
    autoPtr<mapDistribute> mapPtr
    (
        localQueries
        (
            info,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();


    // Do my tests
    // ~~~~~~~~~~~

    const triSurface& s = static_cast<const triSurface&>(*this);

    normal.setSize(triangleIndex.size());

    forAll(triangleIndex, i)
    {
        label trii = triangleIndex[i];
        normal[i] = s[trii].unitNormal(s.points());
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    map.reverseDistribute(info.size(), normal);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getNormal :"
            << " surface " << searchableSurface::name()
            << " finished getting normal for "
            << info.size() << " triangles" << endl;
    }
}


//void Foam::distributedTriSurfaceMesh::getVolumeTypeUncached
//(
//    const pointField& samples,
//    List<volumeType>& volType
//) const
//{
//    if (!Pstream::parRun())
//    {
//        triSurfaceMesh::getVolumeType(samples, volType);
//        return;
//    }
//
//
//    if (!hasVolumeType())
//    {
//        FatalErrorInFunction
//            << "Volume type only supported for closed distributed surfaces."
//            << exit(FatalError);
//    }
//
//    // Trigger (so parallel synchronised) construction of outside type.
//    // Normally this would get triggered from inside individual searches
//    // so would not be parallel synchronised
//    if (outsideVolType_ == volumeType::UNKNOWN)
//    {
//        // Determine nearest (in parallel)
//        const point outsidePt(bounds().max() + 0.5*bounds().span());
//        if (debug)
//        {
//            Pout<< "distributedTriSurfaceMesh::outsideVolumeType :"
//                << " triggering outsidePoint" << outsidePt
//                << " orientation" << endl;
//        }
//
//        const pointField outsidePts(1, outsidePt);
//        List<pointIndexHit> nearestInfo;
//        findNearest
//        (
//            outsidePts,
//            scalarField(1, Foam::sqr(GREAT)),
//            nearestInfo
//        );
//
//        List<volumeType> outsideVolTypes;
//        surfaceSide(outsidePts, nearestInfo, outsideVolTypes);
//        outsideVolType_ = outsideVolTypes[0];
//
//        if (debug)
//        {
//            Pout<< "distributedTriSurfaceMesh::outsideVolumeType :"
//                << " determined outsidePoint" << outsidePt
//                << " to be " << volumeType::names[outsideVolType_] << endl;
//        }
//    }
//
//    // Determine nearest (in parallel)
//    List<pointIndexHit> nearestInfo(samples.size());
//    findNearest
//    (
//        samples,
//        scalarField(samples.size(), Foam::sqr(GREAT)),
//        nearestInfo
//    );
//
//    // Determine orientation (in parallel)
//    surfaceSide(samples, nearestInfo, volType);
//}


void Foam::distributedTriSurfaceMesh::getVolumeType
(
    const pointField& samples,
    List<volumeType>& volType
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::getVolumeType(samples, volType);
        return;
    }


    if (!hasVolumeType())
    {
        FatalErrorInFunction
            << "Volume type only supported for closed distributed surfaces."
            << exit(FatalError);
    }

    // Trigger (so parallel synchronised) construction of outside type.
    // Normally this would get triggered from inside individual searches
    // so would not be parallel synchronised
    if (outsideVolType_ == volumeType::UNKNOWN)
    {
        addProfiling
        (
            getVolumeType,
            "distributedTriSurfaceMesh::getCachedVolumeType"
        );

        // Determine nearest (in parallel)
        const point outsidePt(bounds().max() + 0.5*bounds().span());
        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::getVolumeType :"
                << " surface " << searchableSurface::name()
                << " triggering outsidePoint" << outsidePt
                << " orientation" << endl;
        }

        const pointField outsidePts(1, outsidePt);
        List<pointIndexHit> nearestInfo;
        findNearest
        (
            outsidePts,
            scalarField(1, Foam::sqr(GREAT)),
            nearestInfo
        );

        List<volumeType> outsideVolTypes;
        surfaceSide(outsidePts, nearestInfo, outsideVolTypes);

        // All processors (that have enough surface) will return the same
        // status, all others will return UNKNOWN. Make INSIDE/OUTSIDE win.
        outsideVolType_ = volumeType
        (
            returnReduce(int(outsideVolTypes[0]), maxOp<int>())
        );

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::getVolumeType :"
                << " surface " << searchableSurface::name()
                << " determined outsidePoint" << outsidePt
                << " to be " << volumeType::names[outsideVolType_] << endl;
        }

        if
        (
            outsideVolType_ == volumeType::INSIDE
         || outsideVolType_ == volumeType::OUTSIDE
        )
        {
            // Get local tree
            const indexedOctree<treeDataTriSurface>& t = tree();
            const auto& nodes = t.nodes();
            PackedList<2>& nt = t.nodeTypes();
            nt.resize(nodes.size());
            nt = volumeType::UNKNOWN;

            // Collect midpoints
            DynamicField<point> midPoints(label(0.5*nodes.size()));
            if (nodes.size())
            {
                collectLeafMids(0, midPoints);
            }

            if (debug)
            {
                Pout<< "distributedTriSurfaceMesh::getVolumeType :"
                    << " surface " << searchableSurface::name()
                    << " triggering orientation caching for "
                    << midPoints.size() << " leaf mids" << endl;
            }

            // Get volume type of mid points
            List<volumeType> midVolTypes;

            // Note: could recurse here (since now outsideVolType_ is set)
            // but this would use the cached mid point data which we're trying
            // to calculate. Instead bypass caching and do a full search
            {
                List<pointIndexHit> nearestInfo;
                findNearest
                (
                    midPoints,
                    scalarField(midPoints.size(), Foam::sqr(GREAT)),
                    nearestInfo
                );
                surfaceSide(midPoints, nearestInfo, midVolTypes);
            }

            // Cache on local tree
            if (nodes.size())
            {
                label index = 0;
                calcVolumeType
                (
                    midVolTypes,
                    index,
                    nt,
                    0               // nodeI
                );
            }
            if (debug)
            {
                Pout<< "distributedTriSurfaceMesh::getVolumeType :"
                    << " surface " << searchableSurface::name()
                    << " done orientation caching for "
                    << midPoints.size() << " leaf mids" << endl;
            }
        }
    }


    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getVolumeType :"
            << " surface " << searchableSurface::name()
            << " finding orientation for " << samples.size()
            << " samples" << endl;
    }

    addProfiling
    (
        getVolumeType,
        "distributedTriSurfaceMesh::getVolumeType"
    );


    DynamicList<label> outsideSamples;

    // Distribute samples to relevant processors
    autoPtr<mapDistribute> mapPtr;
    {
        labelListList sendMap(Pstream::nProcs());
        {
            // 1. Count
            labelList nSend(Pstream::nProcs(), Zero);
            forAll(samples, samplei)
            {
                // Find the processors this sample overlaps.
                label nOverlap = 0;
                forAll(procBb_, proci)
                {
                    if (contains(procBb_[proci], samples[samplei]))
                    {
                        nSend[proci]++;
                        nOverlap++;
                    }
                }

                // Special case: point is outside all bbs. These would not
                // get sent to anyone so handle locally. Note that is the
                // equivalent of the test in triSurfaceMesh against the local
                // tree bb
                if (nOverlap == 0)
                {
                    outsideSamples.append(samplei);
                }
            }

            forAll(nSend, proci)
            {
                sendMap[proci].resize_nocopy(nSend[proci]);
                nSend[proci] = 0;
            }

            // 2. Fill
            forAll(samples, samplei)
            {
                // Find the processors this sample overlaps.
                forAll(procBb_, proci)
                {
                    if (contains(procBb_[proci], samples[samplei]))
                    {
                        sendMap[proci][nSend[proci]++] = samplei;
                    }
                }
            }
        }

        mapPtr.reset(new mapDistribute(std::move(sendMap)));
    }
    const mapDistribute& map = *mapPtr;

    // Get the points I need to test
    pointField localPoints(samples);
    map.distribute(localPoints);

    volType.setSize(localPoints.size());
    volType = volumeType::UNKNOWN;

    // Split the local queries into those that I can look up on the tree and
    // those I need to search the nearest for
    DynamicField<point> fullSearchPoints(localPoints.size());
    DynamicList<label> fullSearchMap(localPoints.size());
    forAll(localPoints, i)
    {
        if (tree().nodes().size())
        {
            volType[i] = cachedVolumeType(0, localPoints[i]);
        }
        if (volType[i] == volumeType::UNKNOWN)
        {
            fullSearchMap.append(i);
            fullSearchPoints.append(localPoints[i]);
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getVolumeType :"
            << " surface " << searchableSurface::name()
            << " original samples:" << samples.size()
            << " resulting in local queries:"
            << localPoints.size()
            << " of which cached:" << localPoints.size()-fullSearchPoints.size()
            << endl;
    }

    // Determine nearest (in parallel)
    List<pointIndexHit> nearestInfo;
    findNearest
    (
        fullSearchPoints,
        scalarField(fullSearchPoints.size(), Foam::sqr(GREAT)),
        nearestInfo
    );

    // Determine orientation (in parallel)
    List<volumeType> fullSearchType;
    surfaceSide(fullSearchPoints, nearestInfo, fullSearchType);

    // Combine
    forAll(fullSearchMap, i)
    {
        volType[fullSearchMap[i]] = fullSearchType[i];
    }

    // Send back to originator. In case of multiple answers choose inside or
    // outside
    const volumeType zero(volumeType::UNKNOWN);
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        samples.size(),
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        volType,
        zero,
        volumeCombineOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map.comm()
    );


    //// Combine nearest and volType so we can give nicer error messages
    //List<NearType> nearTypes(volType.size());
    //// For all volType we do have the originating point but not the nearest
    //forAll(localPoints, i)
    //{
    //    nearTypes[i] =
    //        NearType(Pair<point>(localPoints[i], Zero), volType[i]);
    //}
    //// But for all full-search we also have the nearest
    //forAll(fullSearchMap, i)
    //{
    //    nearTypes[fullSearchMap[i]] = NearType
    //    (
    //        Pair<point>(fullSearchPoints[i], nearestInfo[i].hitPoint()),
    //        fullSearchType[i]
    //    );
    //}
    //const NearType zero(Pair<point>(Zero, Zero), volumeType::UNKNOWN);
    //mapDistributeBase::distribute
    //(
    //    Pstream::commsTypes::nonBlocking,
    //    List<labelPair>::null(),
    //    samples.size(),
    //    map.constructMap(),
    //    map.constructHasFlip(),
    //    map.subMap(),
    //    map.subHasFlip(),
    //    nearTypes,
    //    zero,
    //    NearTypeCombineOp(),
    //    noOp(),           // no flipping
    //    UPstream::msgType(),
    //    map.comm()
    //);
    //volType.setSize(nearTypes.size());
    //forAll(nearTypes, i)
    //{
    //    volType[i] = nearTypes[i].second();
    //}

    // Add the points outside the bounding box
    for (label samplei : outsideSamples)
    {
        volType[samplei] = outsideVolType_;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getVolumeType :"
            << " surface " << searchableSurface::name()
            << " finished finding orientation for " << samples.size()
            << " samples" << endl;
    }
}


void Foam::distributedTriSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::getField(info, values);
        return;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getField :"
            << " surface " << searchableSurface::name()
            << " retrieving field for "
            << info.size() << " triangles" << endl;
    }

    addProfiling(getField, "distributedTriSurfaceMesh::getField");

    const auto* fldPtr = findObject<triSurfaceLabelField>("values");

    if (fldPtr)
    {
        const triSurfaceLabelField& fld = *fldPtr;

        // Get query data (= local index of triangle)
        // ~~~~~~~~~~~~~~

        labelList triangleIndex(info.size());
        autoPtr<mapDistribute> mapPtr
        (
            localQueries
            (
                info,
                triangleIndex
            )
        );
        const mapDistribute& map = mapPtr();


        // Do my tests
        // ~~~~~~~~~~~

        values.setSize(triangleIndex.size());

        forAll(triangleIndex, i)
        {
            label trii = triangleIndex[i];
            values[i] = fld[trii];
        }


        // Send back results
        // ~~~~~~~~~~~~~~~~~

        map.reverseDistribute(info.size(), values);
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getField :"
            << " surface " << searchableSurface::name()
            << " finished retrieving field for "
            << info.size() << " triangles" << endl;
    }
}


void Foam::distributedTriSurfaceMesh::overlappingSurface
(
    const triSurface& s,
    const List<treeBoundBox>& bbs,
    boolList& includedFace
)
{
    // Determine what triangles to keep.
    includedFace.setSize(s.size());
    includedFace = false;

    // Create a slightly larger bounding box.
    List<treeBoundBox> bbsX(bbs.size());
    const scalar eps = 1.0e-4;
    forAll(bbs, i)
    {
        const point mid = bbs[i].centre();
        const vector halfSpan = (1.0+eps)*(bbs[i].max() - mid);

        bbsX[i].min() = mid - halfSpan;
        bbsX[i].max() = mid + halfSpan;
    }

    forAll(s, trii)
    {
        const labelledTri& f = s[trii];

        triPointRef tri(s.points(), f);

        if (overlaps(bbsX, tri))
        {
            includedFace[trii] = true;
        }
    }
}


// Subset the part of surface that is overlapping bb.
Foam::triSurface Foam::distributedTriSurfaceMesh::overlappingSurface
(
    const triSurface& s,
    const List<treeBoundBox>& bbs,

    labelList& subPointMap,
    labelList& subFaceMap
)
{
    // Determine what triangles to keep.
    boolList includedFace;
    overlappingSurface(s, bbs, includedFace);
    return subsetMesh(s, includedFace, subPointMap, subFaceMap);
}


// Exchanges indices to the processor they come from.
// - calculates exchange map
// - uses map to calculate local triangle index
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::localQueries
(
    const List<pointIndexHit>& info,
    labelList& triangleIndex
) const
{
    triangleIndex.setSize(info.size());

    const globalIndex& triIndexer = globalTris();


    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // Since determining which processor the query should go to is
    // cheap we do a multi-pass algorithm to save some memory temporarily.

    // 1. Count
    labelList nSend(Pstream::nProcs(), Zero);

    forAll(info, i)
    {
        if (info[i].hit())
        {
            label proci = triIndexer.whichProcID(info[i].index());
            nSend[proci]++;
        }
    }

    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());
    forAll(nSend, proci)
    {
        sendMap[proci].resize_nocopy(nSend[proci]);
        nSend[proci] = 0;
    }

    // 3. Fill sendMap
    forAll(info, i)
    {
        if (info[i].hit())
        {
            label proci = triIndexer.whichProcID(info[i].index());
            triangleIndex[i] = triIndexer.toLocal(proci, info[i].index());
            sendMap[proci][nSend[proci]++] = i;
        }
        else
        {
            triangleIndex[i] = -1;
        }
    }


    // Pack into distribution map
    auto mapPtr = autoPtr<mapDistribute>::New(std::move(sendMap));

    // Send over queries
    mapPtr().distribute(triangleIndex);

    return mapPtr;
}


void Foam::distributedTriSurfaceMesh::distribute
(
    const List<treeBoundBox>& bbs,
    const bool keepNonLocal,
    autoPtr<mapDistribute>& faceMap,
    autoPtr<mapDistribute>& pointMap
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::distribute :"
            << " surface " << searchableSurface::name()
            << " distributing surface according to method:"
            << distributionTypeNames_[distType_]
            << " follow bbs:" << flatOutput(bbs) << endl;
    }

    addProfiling(distribute, "distributedTriSurfaceMesh::distribute");


    // Get bbs of all domains
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        List<List<treeBoundBox>> newProcBb(Pstream::nProcs());

        switch (distType_)
        {
            case FOLLOW:
                newProcBb[Pstream::myProcNo()] = bbs;
                Pstream::allGatherList(newProcBb);
            break;

            case INDEPENDENT:
            case DISTRIBUTED:
                if (currentDistType_ == distType_)
                {
                    return;
                }
                newProcBb = independentlyDistributedBbs(*this);
            break;

            case FROZEN:
                return;
            break;

            default:
                FatalErrorInFunction
                    << "Unsupported distribution type." << exit(FatalError);
            break;
        }

        if (newProcBb == procBb_)
        {
            return;
        }
        else
        {
            procBb_.transfer(newProcBb);
            dict_.set("bounds", procBb_[Pstream::myProcNo()]);
        }
    }


    // Debug information
    if (debug)
    {
        labelList nTris
        (
            UPstream::listGatherValues<label>(triSurface::size())
        );

        if (Pstream::master())
        {
            InfoInFunction
                << "before distribution:" << endl << "\tproc\ttris" << endl;

            forAll(nTris, proci)
            {
                Info<< '\t' << proci << '\t' << nTris[proci] << endl;
            }
            Info<< endl;
        }
    }


    // Use procBbs to determine which faces go where
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList faceSendMap(Pstream::nProcs());
    labelListList pointSendMap(Pstream::nProcs());

    forAll(procBb_, proci)
    {
        overlappingSurface
        (
            *this,
            procBb_[proci],
            pointSendMap[proci],
            faceSendMap[proci]
        );
    }

    if (keepNonLocal)
    {
        // Include in faceSendMap/pointSendMap the triangles that are
        // not mapped to any processor so they stay local.

        const triSurface& s = static_cast<const triSurface&>(*this);

        boolList includedFace(s.size(), true);

        forAll(faceSendMap, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                forAll(faceSendMap[proci], i)
                {
                    includedFace[faceSendMap[proci][i]] = false;
                }
            }
        }

        // Combine includedFace (all triangles that are not on any neighbour)
        // with overlap.

        forAll(faceSendMap[Pstream::myProcNo()], i)
        {
            includedFace[faceSendMap[Pstream::myProcNo()][i]] = true;
        }

        subsetMesh
        (
            s,
            includedFace,
            pointSendMap[Pstream::myProcNo()],
            faceSendMap[Pstream::myProcNo()]
        );
    }


    // Exchange surfaces
    // ~~~~~~~~~~~~~~~~~

    // Storage for resulting surface
    List<labelledTri> allTris;
    pointField allPoints;

    labelListList faceConstructMap(Pstream::nProcs());
    labelListList pointConstructMap(Pstream::nProcs());


    // My own surface first
    // ~~~~~~~~~~~~~~~~~~~~

    {
        labelList pointMap;
        triSurface subSurface
        (
            subsetMesh
            (
                *this,
                faceSendMap[Pstream::myProcNo()],
                pointMap
            )
        );

        allTris = subSurface;
        allPoints = subSurface.points();

        faceConstructMap[Pstream::myProcNo()] = identity
        (
            faceSendMap[Pstream::myProcNo()].size()
        );
        pointConstructMap[Pstream::myProcNo()] = identity
        (
            pointSendMap[Pstream::myProcNo()].size()
        );
    }



    // Send all
    // ~~~~~~~~

    PstreamBuffers pBufs;

    forAll(faceSendMap, proci)
    {
        if (proci != Pstream::myProcNo() && !faceSendMap[proci].empty())
        {
            UOPstream os(proci, pBufs);

            labelList pointMap;
            triSurface subSurface
            (
                subsetMesh
                (
                    *this,
                    faceSendMap[proci],
                    pointMap
                )
            );
            os << subSurface;
        }
    }

    pBufs.finishedSends();


    // Receive and merge all
    // ~~~~~~~~~~~~~~~~~~~~~

    for (const int proci : pBufs.allProcs())
    {
        if (proci != Pstream::myProcNo() && pBufs.recvDataCount(proci))
        {
            UIPstream is(proci, pBufs);

            // Receive
            triSurface subSurface(is);

            // Merge into allSurf
            merge
            (
                mergeDist_,
                subSurface,
                subSurface.points(),

                allTris,
                allPoints,
                faceConstructMap[proci],
                pointConstructMap[proci]
            );
        }
    }


    faceMap.reset
    (
        new mapDistribute
        (
            allTris.size(),
            std::move(faceSendMap),
            std::move(faceConstructMap)
        )
    );
    pointMap.reset
    (
        new mapDistribute
        (
            allPoints.size(),
            std::move(pointSendMap),
            std::move(pointConstructMap)
        )
    );

    // Construct triSurface. Reuse storage.
    triSurface::operator=(triSurface(allTris, patches(), allPoints, true));

    // Clear trees, preserve topological info (surfaceClosed, outsidePointType)
    clearOut();

    // Set the bounds() value to the boundBox of the undecomposed surface
    bounds() = boundBox(points(), true);

    currentDistType_ = distType_;

    // Regions stays same
    // Volume type stays same.

    distributeFields<label>(faceMap());
    distributeFields<scalar>(faceMap());
    distributeFields<vector>(faceMap());
    distributeFields<sphericalTensor>(faceMap());
    distributeFields<symmTensor>(faceMap());
    distributeFields<tensor>(faceMap());

    if (debug)
    {
        labelList nTris
        (
            UPstream::listGatherValues<label>(triSurface::size())
        );

        if (Pstream::master())
        {
            InfoInFunction
                << "after distribution:" << endl << "\tproc\ttris" << endl;

            forAll(nTris, proci)
            {
                Info<< '\t' << proci << '\t' << nTris[proci] << endl;
            }
            Info<< endl;
        }

        if (debug & 2)
        {
            OBJstream str
            (
                searchableSurface::time().path()
               /searchableSurface::name() + "_after.obj"
            );
            Info<< "Writing local bounding box to " << str.name() << endl;
            {
                for (const treeBoundBox& bb : procBb_[Pstream::myProcNo()])
                {
                    str.write(bb, true);  // lines
                }
            }
        }
        if (debug & 2)
        {
            OBJstream str
            (
                searchableSurface::time().path()
               /searchableSurface::name() + "_after_all.obj"
            );
            Info<< "Writing all bounding boxes to " << str.name() << endl;
            for (const auto& myBbs : procBb_)
            {
                for (const treeBoundBox& bb : myBbs)
                {
                    str.write(bb, true);  // lines
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::distribute :"
            << " surface " << searchableSurface::name()
            << " done distributing surface according to method:"
            << distributionTypeNames_[distType_]
            << " follow bbs:" << flatOutput(bbs) << endl;
    }
}


bool Foam::distributedTriSurfaceMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::writeObject :"
            << " surface " << searchableSurface::name()
            << " writing surface:" << writeOnProc << endl;
    }

    // Make sure dictionary goes to same directory as surface
    const_cast<fileName&>(dict_.instance()) = searchableSurface::instance();

    // Copy of triSurfaceMesh::writeObject except for the sorting of
    // triangles by region. This is done so we preserve region names,
    // even if locally we have zero triangles.
    {
        fileName fullPath(searchableSurface::objectPath());

        if (!mkDir(fullPath.path()))
        {
            return false;
        }

        // Important: preserve any zero-sized patches
        triSurface::write(fullPath, true);

        if (!isFile(fullPath))
        {
            return false;
        }
    }

    // Dictionary needs to be written in ascii - binary output not supported.
    streamOpt.format(IOstreamOption::ASCII);
    bool ok = dict_.writeObject(streamOpt, true);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::writeObject :"
            << " surface " << searchableSurface::name()
            << " done writing surface" << endl;
    }

    return ok;
}


void Foam::distributedTriSurfaceMesh::writeStats(Ostream& os) const
{
    boundBox bb;
    label nPoints;
    PatchTools::calcBounds(static_cast<const triSurface&>(*this), bb, nPoints);
    bb.reduce();

    os  << "Triangles    : " << returnReduce(triSurface::size(), sumOp<label>())
        << endl
        << "Vertices     : " << returnReduce(nPoints, sumOp<label>()) << endl
        << "Bounding Box : " << bb << endl
        << "Closed       : " << surfaceClosed_ << endl
        << "Outside point: " << volumeType::names[outsideVolType_] << endl
        << "Distribution : " << distributionTypeNames_[distType_] << endl;
}


// ************************************************************************* //
