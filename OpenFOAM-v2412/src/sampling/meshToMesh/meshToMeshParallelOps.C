/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "meshToMesh.H"
#include "OFstream.H"
#include "Time.H"
#include "globalIndex.H"
#include "mergePoints.H"
#include "processorPolyPatch.H"
#include "SubField.H"
#include "AABBTree.H"
#include "cellBox.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshToMesh::calcDistribution
(
    const polyMesh& src,
    const polyMesh& tgt
) const
{
    label proci = 0;

    if (Pstream::parRun())
    {
        const bitSet hasMesh
        (
            UPstream::listGatherValues<bool>
            (
                src.nCells() > 0 || tgt.nCells() > 0
            )
        );

        const auto nHaveMesh = hasMesh.count();

        if (nHaveMesh == 1)
        {
            proci = hasMesh.find_first();
            DebugInFunction
                << "Meshes local to processor" << proci << endl;
        }
        else if (nHaveMesh > 1)
        {
            proci = -1;
            DebugInFunction
                << "Meshes split across multiple processors" << endl;
        }

        Pstream::broadcast(proci);
    }

    return proci;
}


Foam::label Foam::meshToMesh::calcOverlappingProcs
(
    const List<treeBoundBoxList>& procBb,
    const boundBox& bb,
    boolList& overlaps
) const
{
    overlaps = false;

    label nOverlaps = 0;

    forAll(procBb, proci)
    {
        const treeBoundBoxList& bbp = procBb[proci];

        for (const treeBoundBox& b : bbp)
        {
            if (b.overlaps(bb))
            {
                overlaps[proci] = true;
                ++nOverlaps;
                break;
            }
        }
    }

    return nOverlaps;
}


Foam::autoPtr<Foam::mapDistribute> Foam::meshToMesh::calcProcMap
(
    const polyMesh& src,
    const polyMesh& tgt
) const
{
    // Uses linear construct order
    const
        mapDistributeBase::layoutTypes constructLayout =
        mapDistributeBase::layoutTypes::linear;

    switch (procMapMethod_)
    {
        case procMapMethod::pmLOD:
        {
            Info<< "meshToMesh: Using processorLOD method" << endl;

            // Create processor map of overlapping faces. This map gets
            // (possibly remote) cells from the tgt mesh such that they
            // (together) cover all of the src mesh
            const label nGlobalSrcCells = src.globalData().nTotalCells();
            // Note: minimum content has to be > 1 since otherwise
            //       would keep on splitting. 16 is fairly random choice.
            const label cellsPerBox = max(16, 0.001*nGlobalSrcCells);
            typename processorLODs::cellBox boxLOD
            (
                src.cells(),
                src.faces(),
                src.points(),
                tgt.cells(),
                tgt.faces(),
                tgt.points(),
                cellsPerBox,
                src.nCells()
            );

            return boxLOD.map(constructLayout);
            break;
        }
        default:
        {
            Info<< "meshToMesh: Using AABBTree method" << endl;

            // get decomposition of cells on src mesh
            List<treeBoundBoxList> procBb(Pstream::nProcs());

            if (src.nCells() > 0)
            {
                procBb[Pstream::myProcNo()] = AABBTree<labelList>
                (
                    src.cellPoints(),
                    src.points(),
                    false
                ).boundBoxes();
            }
            else
            {
                procBb[Pstream::myProcNo()] = treeBoundBoxList();
            }

            Pstream::allGatherList(procBb);


            if (debug)
            {
                InfoInFunction
                    << "Determining extent of src mesh per processor:" << nl
                    << "\tproc\tbb" << endl;
                forAll(procBb, proci)
                {
                    Info<< '\t' << proci << '\t' << procBb[proci] << endl;
                }
            }


            // Determine which cells of tgt mesh overlaps src mesh per proc
            labelListList sendMap(Pstream::nProcs());

            {
                // per processor indices into all segments to send
                List<DynamicList<label>> dynSendMap(Pstream::nProcs());
                label iniSize = floor(tgt.nCells()/Pstream::nProcs());

                forAll(dynSendMap, proci)
                {
                    dynSendMap[proci].setCapacity(iniSize);
                }

                // work array - whether src processor bb overlaps the tgt cell
                // bounds
                boolList procBbOverlaps(Pstream::nProcs());

                for (label celli = 0; celli < tgt.nCells(); ++celli)
                {
                    // Bounding box of tgt cell
                    boundBox cellBb(tgt.cellBb(celli));

                    // find the overlapping tgt cells on each src processor
                    (void)calcOverlappingProcs(procBb, cellBb, procBbOverlaps);

                    forAll(procBbOverlaps, proci)
                    {
                        if (procBbOverlaps[proci])
                        {
                            dynSendMap[proci].append(celli);
                        }
                    }
                }

                if (debug)
                {
                    Pout<< "Of my " << tgt.nCells()
                        << " target cells I need to send to:" << nl
                        << "\tproc\tcells" << endl;
                    forAll(dynSendMap, proci)
                    {
                        Pout<< '\t' << proci << '\t'
                            << dynSendMap[proci].size() << endl;
                    }
                }

                // Convert DynamicList -> List
                forAll(sendMap, proci)
                {
                    sendMap[proci].transfer(dynSendMap[proci]);
                }
            }


            return autoPtr<mapDistribute>::New
            (
                constructLayout,
                std::move(sendMap)
            );
            break;
        }
    }
}


void Foam::meshToMesh::distributeCells
(
    const mapDistribute& map,
    const polyMesh& tgtMesh,
    const globalIndex& globalI,
    List<pointField>& points,
    List<label>& nInternalFaces,
    List<faceList>& faces,
    List<labelList>& faceOwner,
    List<labelList>& faceNeighbour,
    List<labelList>& cellIDs,
    List<labelList>& nbrProcIDs,
    List<labelList>& procLocalFaceIDs
) const
{
    points.setSize(Pstream::nProcs());
    nInternalFaces.setSize(Pstream::nProcs(), 0);
    faces.setSize(Pstream::nProcs());
    faceOwner.setSize(Pstream::nProcs());
    faceNeighbour.setSize(Pstream::nProcs());
    cellIDs.setSize(Pstream::nProcs());

    nbrProcIDs.setSize(Pstream::nProcs());;
    procLocalFaceIDs.setSize(Pstream::nProcs());;


    PstreamBuffers pBufs;

    for (const int domain : Pstream::allProcs())
    {
        const labelList& sendElems = map.subMap()[domain];

        if (sendElems.size())
        {
            // reverse cell map
            labelList reverseCellMap(tgtMesh.nCells(), -1);
            forAll(sendElems, subCelli)
            {
                reverseCellMap[sendElems[subCelli]] = subCelli;
            }

            DynamicList<face> subFaces(tgtMesh.nFaces());
            DynamicList<label> subFaceOwner(tgtMesh.nFaces());
            DynamicList<label> subFaceNeighbour(tgtMesh.nFaces());

            DynamicList<label> subNbrProcIDs(tgtMesh.nFaces());
            DynamicList<label> subProcLocalFaceIDs(tgtMesh.nFaces());

            label nInternal = 0;

            // internal faces
            forAll(tgtMesh.faceNeighbour(), facei)
            {
                label own = tgtMesh.faceOwner()[facei];
                label nbr = tgtMesh.faceNeighbour()[facei];
                label subOwn = reverseCellMap[own];
                label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr != -1)
                {
                    nInternal++;

                    if (subOwn < subNbr)
                    {
                        subFaces.append(tgtMesh.faces()[facei]);
                        subFaceOwner.append(subOwn);
                        subFaceNeighbour.append(subNbr);
                        subNbrProcIDs.append(-1);
                        subProcLocalFaceIDs.append(-1);
                    }
                    else
                    {
                        subFaces.append(tgtMesh.faces()[facei].reverseFace());
                        subFaceOwner.append(subNbr);
                        subFaceNeighbour.append(subOwn);
                        subNbrProcIDs.append(-1);
                        subProcLocalFaceIDs.append(-1);
                    }
                }
            }

            // boundary faces for new region
            forAll(tgtMesh.faceNeighbour(), facei)
            {
                label own = tgtMesh.faceOwner()[facei];
                label nbr = tgtMesh.faceNeighbour()[facei];
                label subOwn = reverseCellMap[own];
                label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr == -1)
                {
                    subFaces.append(tgtMesh.faces()[facei]);
                    subFaceOwner.append(subOwn);
                    subFaceNeighbour.append(subNbr);
                    subNbrProcIDs.append(-1);
                    subProcLocalFaceIDs.append(-1);
                }
                else if (subOwn == -1 && subNbr != -1)
                {
                    subFaces.append(tgtMesh.faces()[facei].reverseFace());
                    subFaceOwner.append(subNbr);
                    subFaceNeighbour.append(subOwn);
                    subNbrProcIDs.append(-1);
                    subProcLocalFaceIDs.append(-1);
                }
            }

            // boundary faces of existing region
            forAll(tgtMesh.boundaryMesh(), patchi)
            {
                const polyPatch& pp = tgtMesh.boundaryMesh()[patchi];
                const auto* procPatch = isA<processorPolyPatch>(pp);

                // Store info for faces on processor patches
                const label nbrProci =
                    (procPatch ? procPatch->neighbProcNo() : -1);

                forAll(pp, i)
                {
                    label facei = pp.start() + i;
                    label own = tgtMesh.faceOwner()[facei];

                    if (reverseCellMap[own] != -1)
                    {
                        subFaces.append(tgtMesh.faces()[facei]);
                        subFaceOwner.append(reverseCellMap[own]);
                        subFaceNeighbour.append(-1);
                        subNbrProcIDs.append(nbrProci);
                        subProcLocalFaceIDs.append(i);
                    }
                }
            }

            // reverse point map
            labelList reversePointMap(tgtMesh.nPoints(), -1);
            DynamicList<point> subPoints(tgtMesh.nPoints());
            forAll(subFaces, subFacei)
            {
                face& f = subFaces[subFacei];
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (reversePointMap[pointi] == -1)
                    {
                        reversePointMap[pointi] = subPoints.size();
                        subPoints.append(tgtMesh.points()[pointi]);
                    }

                    f[fp] = reversePointMap[pointi];
                }
            }

            // tgt cells into global numbering
            labelList globalElems(globalI.toGlobal(sendElems));

            if (debug > 1)
            {
                forAll(sendElems, i)
                {
                    Pout<< "tgtProc:" << Pstream::myProcNo()
                        << " sending tgt cell " << sendElems[i]
                        << "[" << globalElems[i] << "]"
                        << " to srcProc " << domain << endl;
                }
            }

            // pass data
            if (domain == Pstream::myProcNo())
            {
                // allocate my own data
                points[Pstream::myProcNo()] = subPoints;
                nInternalFaces[Pstream::myProcNo()] = nInternal;
                faces[Pstream::myProcNo()] = subFaces;
                faceOwner[Pstream::myProcNo()] = subFaceOwner;
                faceNeighbour[Pstream::myProcNo()] = subFaceNeighbour;
                cellIDs[Pstream::myProcNo()] = globalElems;
                nbrProcIDs[Pstream::myProcNo()] = subNbrProcIDs;
                procLocalFaceIDs[Pstream::myProcNo()] = subProcLocalFaceIDs;
            }
            else
            {
                // send data to other processor domains
                UOPstream toDomain(domain, pBufs);

                toDomain
                    << subPoints
                    << nInternal
                    << subFaces
                    << subFaceOwner
                    << subFaceNeighbour
                    << globalElems
                    << subNbrProcIDs
                    << subProcLocalFaceIDs;
            }
        }
    }

    // Start receiving
    pBufs.finishedSends();

    // Consume
    for (const int domain : Pstream::allProcs())
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo() && recvElems.size())
        {
            UIPstream str(domain, pBufs);

            str >> points[domain]
                >> nInternalFaces[domain]
                >> faces[domain]
                >> faceOwner[domain]
                >> faceNeighbour[domain]
                >> cellIDs[domain]
                >> nbrProcIDs[domain]
                >> procLocalFaceIDs[domain];
        }

        if (debug)
        {
            Pout<< "Target mesh send sizes[" << domain << "]"
                << ": points="<< points[domain].size()
                << ", faces=" << faces[domain].size()
                << ", nInternalFaces=" << nInternalFaces[domain]
                << ", faceOwn=" << faceOwner[domain].size()
                << ", faceNbr=" << faceNeighbour[domain].size()
                << ", cellIDs=" << cellIDs[domain].size() << endl;
        }
    }
}


void Foam::meshToMesh::distributeAndMergeCells
(
    const mapDistribute& map,
    const polyMesh& tgt,
    const globalIndex& globalI,
    pointField& tgtPoints,
    faceList& tgtFaces,
    labelList& tgtFaceOwners,
    labelList& tgtFaceNeighbours,
    labelList& tgtCellIDs
) const
{
    // Exchange per-processor data
    List<pointField> allPoints;
    List<label> allNInternalFaces;
    List<faceList> allFaces;
    List<labelList> allFaceOwners;
    List<labelList> allFaceNeighbours;
    List<labelList> allTgtCellIDs;

    // Per target mesh face the neighbouring proc and index in
    // processor patch (all -1 for normal boundary face)
    List<labelList> allNbrProcIDs;
    List<labelList> allProcLocalFaceIDs;

    distributeCells
    (
        map,
        tgt,
        globalI,
        allPoints,
        allNInternalFaces,
        allFaces,
        allFaceOwners,
        allFaceNeighbours,
        allTgtCellIDs,
        allNbrProcIDs,
        allProcLocalFaceIDs
    );

    // Convert lists into format that can be used to generate a valid polyMesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Points and cells are collected into single flat lists:
    // - i.e. proc0, proc1 ... procN
    //
    // Faces need to be sorted after collection to that internal faces are
    // contiguous, followed by all boundary faces
    //
    // Processor patch faces between included cells on neighbouring processors
    // are converted into internal faces
    //
    // Face list structure:
    // - Per processor:
    //   - internal faces
    //   - processor faces that have been converted into internal faces
    // - Followed by all boundary faces
    //   - from 'normal' boundary faces
    //   - from singularly-sided processor patch faces


    // Number of internal+coupled faces
    labelList allNIntCoupledFaces(allNInternalFaces);

    // Starting offset for points
    label nPoints = 0;
    labelList pointOffset(Pstream::nProcs(), Zero);
    forAll(allPoints, proci)
    {
        pointOffset[proci] = nPoints;
        nPoints += allPoints[proci].size();
    }

    // Starting offset for cells
    label nCells = 0;
    labelList cellOffset(Pstream::nProcs(), Zero);
    forAll(allTgtCellIDs, proci)
    {
        cellOffset[proci] = nCells;
        nCells += allTgtCellIDs[proci].size();
    }

    // Count any coupled faces
    typedef FixedList<label, 3> label3;
    typedef HashTable<label, label3> procCoupleInfo;
    procCoupleInfo procFaceToGlobalCell;

    forAll(allNbrProcIDs, proci)
    {
        const labelList& nbrProci = allNbrProcIDs[proci];
        const labelList& localFacei = allProcLocalFaceIDs[proci];

        forAll(nbrProci, i)
        {
            if (nbrProci[i] != -1 && localFacei[i] != -1)
            {
                label3 key;
                key[0] = min(proci, nbrProci[i]);
                key[1] = max(proci, nbrProci[i]);
                key[2] = localFacei[i];

                const auto fnd = procFaceToGlobalCell.cfind(key);

                if (!fnd.good())
                {
                    procFaceToGlobalCell.insert(key, -1);
                }
                else
                {
                    if (debug > 1)
                    {
                        Pout<< "Additional internal face between procs:"
                            << key[0] << " and " << key[1]
                            << " across local face " << key[2] << endl;
                    }

                    allNIntCoupledFaces[proci]++;
                }
            }
        }
    }


    // Starting offset for internal faces
    label nIntFaces = 0;
    label nFacesTotal = 0;
    labelList internalFaceOffset(Pstream::nProcs(), Zero);
    forAll(allNIntCoupledFaces, proci)
    {
        label nCoupledFaces =
            allNIntCoupledFaces[proci] - allNInternalFaces[proci];

        internalFaceOffset[proci] = nIntFaces;
        nIntFaces += allNIntCoupledFaces[proci];
        nFacesTotal += allFaceOwners[proci].size() - nCoupledFaces;
    }

    tgtPoints.setSize(nPoints);
    tgtFaces.setSize(nFacesTotal);
    tgtFaceOwners.setSize(nFacesTotal);
    tgtFaceNeighbours.setSize(nFacesTotal);
    tgtCellIDs.setSize(nCells);

    // Insert points
    forAll(allPoints, proci)
    {
        const pointField& pts = allPoints[proci];
        SubList<point>(tgtPoints, pts.size(), pointOffset[proci]) = pts;
    }

    // Insert cellIDs
    forAll(allTgtCellIDs, proci)
    {
        const labelList& cellIDs = allTgtCellIDs[proci];
        SubList<label>(tgtCellIDs, cellIDs.size(), cellOffset[proci]) = cellIDs;
    }


    // Insert internal faces (from internal faces)
    forAll(allFaces, proci)
    {
        const faceList& fcs = allFaces[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const labelList& faceNs = allFaceNeighbours[proci];

        SubList<face> slice
        (
            tgtFaces,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        slice = SubList<face>(fcs, allNInternalFaces[proci]);
        forAll(slice, i)
        {
            add(slice[i], pointOffset[proci]);
        }

        SubField<label> ownSlice
        (
            tgtFaceOwners,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        ownSlice = SubField<label>(faceOs, allNInternalFaces[proci]);
        add(ownSlice, cellOffset[proci]);

        SubField<label> nbrSlice
        (
            tgtFaceNeighbours,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        nbrSlice = SubField<label>(faceNs, allNInternalFaces[proci]);
        add(nbrSlice, cellOffset[proci]);

        internalFaceOffset[proci] += allNInternalFaces[proci];
    }


    // Insert internal faces (from coupled face-pairs)
    forAll(allNbrProcIDs, proci)
    {
        const labelList& nbrProci = allNbrProcIDs[proci];
        const labelList& localFacei = allProcLocalFaceIDs[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const faceList& fcs = allFaces[proci];

        forAll(nbrProci, i)
        {
            if (nbrProci[i] != -1 && localFacei[i] != -1)
            {
                label3 key;
                key[0] = min(proci, nbrProci[i]);
                key[1] = max(proci, nbrProci[i]);
                key[2] = localFacei[i];

                auto fnd = procFaceToGlobalCell.find(key);

                if (fnd.good())
                {
                    label tgtFacei = fnd();
                    if (tgtFacei == -1)
                    {
                        // on first visit store the new cell on this side
                        fnd() = cellOffset[proci] + faceOs[i];
                    }
                    else
                    {
                        // get owner and neighbour in new cell numbering
                        label newOwn = cellOffset[proci] + faceOs[i];
                        label newNbr = fnd();
                        label tgtFacei = internalFaceOffset[proci]++;

                        if (debug > 1)
                        {
                            Pout<< "    proc " << proci
                                << "\tinserting face:" << tgtFacei
                                << " connection between owner " << newOwn
                                << " and neighbour " << newNbr
                                << endl;
                        }

                        if (newOwn < newNbr)
                        {
                            // we have correct orientation
                            tgtFaces[tgtFacei] = fcs[i];
                            tgtFaceOwners[tgtFacei] = newOwn;
                            tgtFaceNeighbours[tgtFacei] = newNbr;
                        }
                        else
                        {
                            // reverse orientation
                            tgtFaces[tgtFacei] = fcs[i].reverseFace();
                            tgtFaceOwners[tgtFacei] = newNbr;
                            tgtFaceNeighbours[tgtFacei] = newOwn;
                        }

                        add(tgtFaces[tgtFacei], pointOffset[proci]);

                        // mark with unique value
                        fnd() = -2;
                    }
                }
            }
        }
    }


    forAll(allNbrProcIDs, proci)
    {
        const labelList& nbrProci = allNbrProcIDs[proci];
        const labelList& localFacei = allProcLocalFaceIDs[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const labelList& faceNs = allFaceNeighbours[proci];
        const faceList& fcs = allFaces[proci];

        forAll(nbrProci, i)
        {
            // coupled boundary face
            if (nbrProci[i] != -1 && localFacei[i] != -1)
            {
                label3 key;
                key[0] = min(proci, nbrProci[i]);
                key[1] = max(proci, nbrProci[i]);
                key[2] = localFacei[i];

                label tgtFacei = procFaceToGlobalCell[key];

                if (tgtFacei == -1)
                {
                    FatalErrorInFunction
                        << "Unvisited " << key
                        << abort(FatalError);
                }
                else if (tgtFacei != -2)
                {
                    label newOwn = cellOffset[proci] + faceOs[i];
                    label tgtFacei = nIntFaces++;

                    if (debug > 1)
                    {
                        Pout<< "    proc " << proci
                            << "\tinserting boundary face:" << tgtFacei
                            << " from coupled face " << key
                            << endl;
                    }

                    tgtFaces[tgtFacei] = fcs[i];
                    add(tgtFaces[tgtFacei], pointOffset[proci]);

                    tgtFaceOwners[tgtFacei] = newOwn;
                    tgtFaceNeighbours[tgtFacei] = -1;
                }
            }
            // normal boundary face
            else
            {
                label own = faceOs[i];
                label nbr = faceNs[i];
                if ((own != -1) && (nbr == -1))
                {
                    label newOwn = cellOffset[proci] + faceOs[i];
                    label tgtFacei = nIntFaces++;

                    tgtFaces[tgtFacei] = fcs[i];
                    add(tgtFaces[tgtFacei], pointOffset[proci]);

                    tgtFaceOwners[tgtFacei] = newOwn;
                    tgtFaceNeighbours[tgtFacei] = -1;
                }
            }
        }
    }


    if (debug)
    {
        // only merging points in debug mode

        labelList oldToNew;
        label nChanged = Foam::inplaceMergePoints
        (
            tgtPoints,
            SMALL,
            false,
            oldToNew
        );

        if (nChanged)
        {
            Pout<< "Merged from " << oldToNew.size()
                << " down to " << tgtPoints.size() << " points" << endl;

            for (auto& f : tgtFaces)
            {
                inplaceRenumber(oldToNew, f);
            }
        }
    }
}


// ************************************************************************* //
