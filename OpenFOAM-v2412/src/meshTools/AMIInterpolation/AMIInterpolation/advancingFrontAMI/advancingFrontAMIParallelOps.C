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

#include "advancingFrontAMI.H"
#include "mergePoints.H"
#include "mapDistribute.H"
#include "AABBTree.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::advancingFrontAMI::calcOverlappingProcs
(
    const List<treeBoundBoxList>& procBb,
    const treeBoundBox& bb,
    boolList& overlaps
) const
{
    overlaps.setSize(procBb.size());
    overlaps = false;

    label nOverlaps = 0;

    forAll(procBb, proci)
    {
        const treeBoundBoxList& bbp = procBb[proci];

        for (const treeBoundBox& tbb: bbp)
        {
            if (tbb.overlaps(bb))
            {
                overlaps[proci] = true;
                ++nOverlaps;
                break;
            }
        }
    }

    return nOverlaps;
}


void Foam::advancingFrontAMI::distributePatches
(
    const mapDistribute& map,
    const primitivePatch& pp,
    const globalIndex& gi,
    List<faceList>& faces,
    List<pointField>& points,
    List<labelList>& faceIDs
) const
{
    faces.setSize(Pstream::nProcs(comm_));
    points.setSize(Pstream::nProcs(comm_));
    faceIDs.setSize(Pstream::nProcs(comm_));

    PstreamBuffers pBufs(comm_);

    for (const int domain : Pstream::allProcs(comm_))
    {
        const labelList& sendElems = map.subMap()[domain];

        if (sendElems.empty())
        {
            // Safety
            faces[domain].clear();
            points[domain].clear();
            faceIDs[domain].clear();
        }
        else
        {
            faceList subFaces(UIndirectList<face>(pp, sendElems));
            primitivePatch subPatch(SubList<face>(subFaces), pp.points());

            if (debug & 2)
            {
                Pout<< "distributePatches: to processor " << domain
                    << " sending faces " << subPatch.faceCentres() << endl;
            }


            if (domain == Pstream::myProcNo(comm_))
            {
                // Do send/receive for myself
                faces[domain] = subPatch.localFaces();
                points[domain] = subPatch.localPoints();
                faceIDs[domain] =
                    gi.toGlobal(Pstream::myProcNo(comm_), sendElems);
            }
            else
            {
                // Normal send
                UOPstream str(domain, pBufs);
                str
                    << subPatch.localFaces()
                    << subPatch.localPoints()
                    << gi.toGlobal(Pstream::myProcNo(comm_), sendElems);
            }
        }
    }

    pBufs.finishedSends();


    // Consume
    for (const int domain : Pstream::allProcs(comm_))
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo(comm_) && recvElems.size())
        {
            UIPstream is(domain, pBufs);

            is  >> faces[domain]
                >> points[domain]
                >> faceIDs[domain];
        }
    }
}


void Foam::advancingFrontAMI::distributeAndMergePatches
(
    const mapDistribute& map,
    const primitivePatch& tgtPatch,
    const globalIndex& gi,
    faceList& tgtFaces,
    pointField& tgtPoints,
    labelList& tgtFaceIDs
) const
{
    // Exchange per-processor data
    List<faceList> allFaces;
    List<pointField> allPoints;
    List<labelList> allTgtFaceIDs;
    distributePatches(map, tgtPatch, gi, allFaces, allPoints, allTgtFaceIDs);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, proci)
    {
        nFaces += allFaces[proci].size();
        nPoints += allPoints[proci].size();
    }

    tgtFaces.setSize(nFaces);
    tgtPoints.setSize(nPoints);
    tgtFaceIDs.setSize(nFaces);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const labelList& faceIDs = allTgtFaceIDs[Pstream::myProcNo(comm_)];
        SubList<label>(tgtFaceIDs, faceIDs.size()) = faceIDs;

        const faceList& fcs = allFaces[Pstream::myProcNo(comm_)];
        for (const face& f : fcs)
        {
            face& newF = tgtFaces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo(comm_)];
        for (const point& pt: pts)
        {
            tgtPoints[nPoints++] = pt;
        }
    }


    // Other proc data follows
    forAll(allFaces, proci)
    {
        if (proci != Pstream::myProcNo(comm_))
        {
            const labelList& faceIDs = allTgtFaceIDs[proci];
            SubList<label>(tgtFaceIDs, faceIDs.size(), nFaces) = faceIDs;

            const faceList& fcs = allFaces[proci];
            for (const face& f : fcs)
            {
                face& newF = tgtFaces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[proci];
            for (const point& pt: pts)
            {
                tgtPoints[nPoints++] = pt;
            }
        }
    }

    // Merge
    labelList oldToNew;
    bool nChanged = Foam::inplaceMergePoints
    (
        tgtPoints,
        SMALL,
        false,
        oldToNew
    );

    if (nChanged)
    {
        if (debug)
        {
            Pout<< "Merged from " << oldToNew.size()
                << " down to " << tgtPoints.size() << " points" << endl;
        }

        for (face& f : tgtFaces)
        {
            inplaceRenumber(oldToNew, f);
        }
    }
}


Foam::autoPtr<Foam::mapDistribute> Foam::advancingFrontAMI::calcProcMap
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
) const
{
    // Get decomposition of patch
    List<treeBoundBoxList> procBb(Pstream::nProcs(comm_));

    if (srcPatch.size())
    {
        procBb[Pstream::myProcNo(comm_)] =
            AABBTree<face>
            (
                srcPatch.localFaces(),
                srcPatch.localPoints(),
                false
            ).boundBoxes();
    }
    else
    {
        procBb[Pstream::myProcNo(comm_)] = treeBoundBoxList();
    }

    Pstream::allGatherList(procBb, UPstream::msgType(), comm_);

    if (debug)
    {
        Info<< "Determining extent of srcPatch per processor:" << nl
            << "\tproc\tbb" << endl;
        forAll(procBb, proci)
        {
            Info<< '\t' << proci << '\t' << procBb[proci] << endl;
        }
    }

    // Determine which faces of tgtPatch overlaps srcPatch per proc
    const faceList& faces = tgtPatch.localFaces();
    const pointField& points = tgtPatch.localPoints();

    labelListList sendMap;

    {
        // Per processor indices into all segments to send
        List<DynamicList<label>> dynSendMap(Pstream::nProcs(comm_));

        // Work array - whether processor bb overlaps the face bounds
        boolList procBbOverlaps(Pstream::nProcs(comm_));

        forAll(faces, facei)
        {
            if (faces[facei].size())
            {
                treeBoundBox faceBb(points, faces[facei]);

                // Find the processor this face overlaps
                calcOverlappingProcs(procBb, faceBb, procBbOverlaps);

                forAll(procBbOverlaps, proci)
                {
                    if (procBbOverlaps[proci])
                    {
                        dynSendMap[proci].append(facei);
                    }
                }
            }
        }

        // Convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs(comm_));
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }
    }

    // Debug printing
    if (debug)
    {
        Pout<< "Of my " << faces.size() << " I need to send to:" << nl
            << "\tproc\tfaces" << endl;
        forAll(sendMap, proci)
        {
            Pout<< '\t' << proci << '\t' << sendMap[proci].size() << endl;
        }
    }

    return autoPtr<mapDistribute>::New
    (
        std::move(sendMap),
        false,      //subHasFlip
        false,      //constructHasFlip
        comm_
    );
}


// ************************************************************************* //
