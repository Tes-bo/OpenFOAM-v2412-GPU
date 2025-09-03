/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022,2024 OpenCFD Ltd.
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

#include "addPatchCellLayer.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyAddCell.H"
#include "globalIndex.H"
#include "PatchTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(addPatchCellLayer, 0);

    // Reduction class to get maximum value over coupled face.
    class combineEqOpFace
    {
    public:

        void operator()(face& x, const face& y) const
        {
            if (x.size())
            {
                if (y.size())
                {
                    if (x.size() != y.size())
                    {
                        FatalErrorInFunction
                            << "face x:" << flatOutput(x)
                            << " face y:" << flatOutput(y)
                            << exit(FatalError);
                    }

                    label j = 0;
                    forAll(x, i)
                    {
                        x[i] = max(x[i], y[j]);

                        j = y.rcIndex(j);
                    }
                }
            }
            else if (y.size())
            {
                x.setSize(y.size());
                label j = 0;
                forAll(x, i)
                {
                    x[i] = y[j];
                    j = y.rcIndex(j);
                }
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::addPatchCellLayer::nbrFace
(
    const labelListList& edgeFaces,
    const label edgei,
    const label facei
)
{
    const labelList& eFaces = edgeFaces[edgei];

    if (eFaces.size() == 2)
    {
        return (eFaces[0] != facei ? eFaces[0] : eFaces[1]);
    }
    else
    {
        return -1;
    }
}


void Foam::addPatchCellLayer::addVertex
(
    const label pointi,
    face& f,
    label& fp
)
{
    if (fp == 0)
    {
        f[fp++] = pointi;
    }
    else
    {
        if (f[fp-1] != pointi && f[0] != pointi)
        {
            f[fp++] = pointi;
        }
    }
}


// Is edge to the same neighbour? (and needs extrusion and has not been
// dealt with already)
bool Foam::addPatchCellLayer::sameEdgeNeighbour
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const boolList& doneEdge,
    const label thisGlobalFacei,
    const label nbrGlobalFacei,
    const label edgei
) const
{
    const edge& e = pp.edges()[edgei];

    return
        !doneEdge[edgei]                            // not yet handled
     && (
            addedPoints_[e[0]].size()               // is extruded
         || addedPoints_[e[1]].size()
        )
     && (
            nbrFace(globalEdgeFaces, edgei, thisGlobalFacei)
         == nbrGlobalFacei  // is to same neighbour
        );
}


// Collect consecutive string of edges that connects the same two
// (possibly coupled) faces. Returns -1 if no unvisited edge can be found.
// Otherwise returns start and end index in face.
Foam::labelPair Foam::addPatchCellLayer::getEdgeString
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const boolList& doneEdge,
    const label patchFacei,
    const label globalFacei
) const
{
    const labelList& fEdges = pp.faceEdges()[patchFacei];

    label startFp = -1;
    label endFp = -1;

    // Get edge that hasn't been done yet but needs extrusion
    forAll(fEdges, fp)
    {
        label edgei = fEdges[fp];
        const edge& e = pp.edges()[edgei];

        if
        (
            !doneEdge[edgei]
         && ( addedPoints_[e[0]].size() || addedPoints_[e[1]].size() )
        )
        {
            startFp = fp;
            break;
        }
    }

    if (startFp != -1)
    {
        // We found an edge that needs extruding but hasn't been done yet.
        // Now find the face on the other side
        label nbrGlobalFacei = nbrFace
        (
            globalEdgeFaces,
            fEdges[startFp],
            globalFacei
        );

        if (nbrGlobalFacei == -1)
        {
            // Proper boundary edge. Only extrude single edge.
            endFp = startFp;
        }
        else
        {
            // Search back for edge
            // - which hasn't been handled yet
            // - with same neighbour
            // - that needs extrusion

            const label initFp = startFp;
            while (true)
            {
                label prevFp = fEdges.rcIndex(startFp);

                if (prevFp == initFp)
                {
                    const edge& e = pp.edges()[fEdges[initFp]];
                    const face& localF = pp.localFaces()[patchFacei];

                    FatalErrorInFunction
                        << "On face:" << patchFacei
                        << " fc:" << pp.faceCentres()[patchFacei]
                        << " vertices:" << localF
                        << " points:"
                        << UIndirectList<point>(pp.points(), pp[patchFacei])
                        << " edges:" << fEdges
                        << " All edges of face seem to have same neighbour "
                        << nbrGlobalFacei
                        << " starting walking from edge " << e
                        << exit(FatalError);
                }

                if
                (
                    !sameEdgeNeighbour
                    (
                        pp,
                        globalEdgeFaces,
                        doneEdge,
                        globalFacei,
                        nbrGlobalFacei,
                        fEdges[prevFp]
                    )
                )
                {
                    break;
                }
                startFp = prevFp;
            }

            // Search forward for end of string
            endFp = startFp;
            while (true)
            {
                label nextFp = fEdges.fcIndex(endFp);

                if
                (
                    !sameEdgeNeighbour
                    (
                        pp,
                        globalEdgeFaces,
                        doneEdge,
                        globalFacei,
                        nbrGlobalFacei,
                        fEdges[nextFp]
                    )
                )
                {
                    break;
                }
                endFp = nextFp;
            }
        }
    }

    return labelPair(startFp, endFp);
}


Foam::label Foam::addPatchCellLayer::addSideFace
(
    const indirectPrimitivePatch& pp,
    const labelListList& addedCells,    // per pp face the new extruded cell
    const face& newFace,
    const label newPatchID,
    const label zonei,
    const bool edgeFlip,
    const label inflateFacei,

    const label ownFacei,               // pp face that provides owner
    const label nbrFacei,
    const label meshEdgei,              // corresponding mesh edge
    const label layeri,                 // layer
    const label numEdgeFaces,           // number of layers for edge
    const labelList& meshFaces,         // precalculated edgeFaces
    polyTopoChange& meshMod
) const
{
    // Adds a side face i.e. extrudes a patch edge.

    label addedFacei = -1;


    // Is patch edge external edge of indirectPrimitivePatch?
    if (nbrFacei == -1)
    {
        // External edge so external face.

        // Determine if different number of layer on owner and neighbour side
        // (relevant only for coupled faces). See section for internal edge
        // below.

        label layerOwn;

        if (addedCells[ownFacei].size() < numEdgeFaces)
        {
            label offset = numEdgeFaces - addedCells[ownFacei].size();
            if (layeri <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layeri - offset;
            }
        }
        else
        {
            layerOwn = layeri;
        }


        //Pout<< "Added boundary face:" << newFace
        //    << " atfc:" << newFace.centre(meshMod.points())
        //    << " n:" << newFace.unitNormal(meshMod.points())
        //    << " own:" << addedCells[ownFacei][layerOwn]
        //    << " patch:" << newPatchID
        //    << endl;

        addedFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                        // face
                addedCells[ownFacei][layerOwn], // owner
                -1,                             // neighbour
                -1,                             // master point
                -1,                             // master edge
                inflateFacei,                   // master face
                false,                          // flux flip
                newPatchID,                     // patch for face
                zonei,                          // zone for face
                edgeFlip                        // face zone flip
            )
        );
    }
    else
    {
        // When adding side faces we need to modify neighbour and owners
        // in region where layer mesh is stopped. Determine which side
        // has max number of faces and make sure layers match closest to
        // original pp if there are different number of layers.

        label layerNbr;
        label layerOwn;

        if (addedCells[ownFacei].size() > addedCells[nbrFacei].size())
        {
            label offset =
                addedCells[ownFacei].size() - addedCells[nbrFacei].size();

            layerOwn = layeri;

            if (layeri <= offset)
            {
                layerNbr = 0;
            }
            else
            {
                layerNbr = layeri - offset;
            }
        }
        else if (addedCells[nbrFacei].size() > addedCells[ownFacei].size())
        {
            label offset =
                addedCells[nbrFacei].size() - addedCells[ownFacei].size();

            layerNbr = layeri;

            if (layeri <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layeri - offset;
            }
        }
        else
        {
            // Same number of layers on both sides.
            layerNbr = layeri;
            layerOwn = layeri;
        }


        // Check mesh internal faces using edge to initialise
        label inflateEdgei = -1;
        if (addToMesh_)
        {
            forAll(meshFaces, i)
            {
                if (mesh_.isInternalFace(meshFaces[i]))
                {
                    // meshEdge uses internal faces so ok to inflate from it
                    inflateEdgei = meshEdgei;
                    break;
                }
            }
        }


        addedFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                        // face
                addedCells[ownFacei][layerOwn], // owner
                addedCells[nbrFacei][layerNbr], // neighbour
                -1,                             // master point
                inflateEdgei,                   // master edge
                -1,                             // master face
                false,                          // flux flip
                -1,                             // patch for face
                zonei,                          // zone for face
                edgeFlip                        // face zone flip
            )
        );

        //Pout<< "Added internal face:" << newFace
        //    << " atfc:" << newFace.centre(meshMod.points())
        //    << " n:" << newFace.unitNormal(meshMod.points())
        //    << " own:" << addedCells[ownFacei][layerOwn]
        //    << " nei:" << addedCells[nbrFacei][layerNbr]
        //    << endl;
    }

    return addedFacei;
}


void Foam::addPatchCellLayer::setFaceProps
(
    const polyMesh& mesh,
    const label facei,

    label& patchi,
    label& zonei,
    bool& zoneFlip
)
{
    patchi = mesh.boundaryMesh().whichPatch(facei);
    zonei = mesh.faceZones().whichZone(facei);
    zoneFlip = false;
    if (zonei != -1)
    {
        label index = mesh.faceZones()[zonei].whichFace(facei);
        zoneFlip = mesh.faceZones()[zonei].flipMap()[index];
    }
}


void Foam::addPatchCellLayer::setFaceProps
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const label ppEdgeI,
    const label faceI,

    label& patchI,
    label& zoneI,
    bool& zoneFlip,
    label& inflateFaceI
)
{
    setFaceProps
    (
        mesh,
        faceI,

        patchI,
        zoneI,
        zoneFlip
    );

    if (patchI != -1 || zoneI != -1)
    {
        inflateFaceI = faceI;
    }

    if (zoneI != -1)
    {
        // Correct flip for patch edge ordering
        const edge& ppEdge = pp.edges()[ppEdgeI];
        const edge patchEdge
        (
            pp.meshPoints()[ppEdge[0]],
            pp.meshPoints()[ppEdge[1]]
        );

        const face& f = mesh.faces()[faceI];
        bool found = false;
        forAll(f, fp)
        {
            const edge e(f[fp], f.nextLabel(fp));
            int stat = edge::compare(e, patchEdge);
            if (stat == 1)
            {
                found = true;
                break;
            }
            else if (stat == -1)
            {
                found = true;
                zoneFlip = !zoneFlip;
                break;
            }
        }

        if (!found)
        {
            //FatalErrorInFunction
            WarningInFunction
                << "Problem: cannot find patch edge " << ppEdgeI
                << " with mesh vertices " << patchEdge
                << " at " << patchEdge.line(mesh.points())
                << " in face " << faceI << " with mesh vertices "
                << f
                << " at " << pointField(mesh.points(), f)
                << endl
                << "Continuing with potentially incorrect faceZone orientation"
                //<< exit(FatalError);
                << endl;
        }
    }
}


//void Foam::addPatchCellLayer::findZoneFace
//(
//    const bool useInternalFaces,
//    const bool useBoundaryFaces,
//
//    const polyMesh& mesh,
//    const indirectPrimitivePatch& pp,
//    const label ppEdgeI,
//    const labelUIndList& excludeFaces,
//    const labelList& meshFaces,
//
//    label& inflateFaceI,
//    label& patchI,
//    label& zoneI,
//    bool& zoneFlip
//)
//{
//    inflateFaceI = -1;
//    patchI = -1;
//    zoneI = -1;
//    zoneFlip = false;
//
//    forAll(meshFaces, k)
//    {
//        label faceI = meshFaces[k];
//
//        if
//        (
//            !excludeFaces.found(faceI)
//         && (
//                (mesh.isInternalFace(faceI) && useInternalFaces)
//             || (!mesh.isInternalFace(faceI) && useBoundaryFaces)
//            )
//        )
//        {
//            setFaceProps
//            (
//                mesh,
//                pp,
//                ppEdgeI,
//                faceI,
//
//                patchI,
//                zoneI,
//                zoneFlip,
//                inflateFaceI
//            );
//
//            if (zoneI != -1 || patchI != -1)
//            {
//                break;
//            }
//        }
//    }
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::addPatchCellLayer::addPatchCellLayer
(
    const polyMesh& mesh,
    const bool addToMesh,
    const bool extrude
)
:
    mesh_(mesh),
    addToMesh_(addToMesh),
    extrude_(extrude),
    addedPoints_(0),
    layerFaces_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::addPatchCellLayer::addedCells
(
    const polyMesh& mesh,
    const labelListList& layerFaces
)
{
    labelListList layerCells(layerFaces.size());

    forAll(layerFaces, patchFacei)
    {
        const labelList& faceLabels = layerFaces[patchFacei];

        if (faceLabels.size())
        {
            labelList& added = layerCells[patchFacei];
            added.setSize(faceLabels.size()-1);

            for (label i = 0; i < faceLabels.size()-1; i++)
            {
                added[i] = mesh.faceOwner()[faceLabels[i+1]];
            }
        }
    }
    return layerCells;
}


Foam::labelListList Foam::addPatchCellLayer::addedCells() const
{
    return addedCells(mesh_, layerFaces_);
}


Foam::labelListList Foam::addPatchCellLayer::globalEdgeFaces
(
    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const indirectPrimitivePatch& pp,
    const bitSet& ppFlip
)
{

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    // From mesh edge to global face labels. Non-empty sublists only for
    // pp edges.
    labelListList globalEdgeFaces(mesh.nEdges());

    const labelListList& edgeFaces = pp.edgeFaces();

    DynamicList<label> faceIDs;

    forAll(edgeFaces, edgeI)
    {
        // Store face and processor as unique tag.
        faceIDs.clear();
        for (const label patchFacei : edgeFaces[edgeI])
        {
            const label facei = pp.addressing()[patchFacei];

            // Ignore if ppFlip is on the other side of the (coupled) boundary
            if
            (
                ppFlip.empty()
            || !ppFlip[patchFacei]
            || mesh.isInternalFace(facei)
            )
            {
                faceIDs.append(globalFaces.toGlobal(facei));
            }
        }

        globalEdgeFaces[meshEdges[edgeI]] = std::move(faceIDs);
    }


    // Synchronise across coupled edges.
    syncTools::syncEdgeList
    (
        mesh,
        globalEdgeFaces,
        ListOps::uniqueEqOp<label>(),
        labelList()             // null value
    );

    // Extract pp part
    return labelListList(UIndirectList<labelList>(globalEdgeFaces, meshEdges));
}


Foam::labelListList Foam::addPatchCellLayer::globalEdgeFaces
(
    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const indirectPrimitivePatch& pp
)
{
    return globalEdgeFaces(mesh, globalFaces, pp, bitSet::null());
}


void Foam::addPatchCellLayer::markPatchEdges
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const labelListList& edgeGlobalFaces,
    const labelList& meshEdges,

    bitSet& isPatchEdge,
    bitSet& isPatchBoundaryEdge
)
{
    // Mark (mesh) edges:
    // - anywhere on extrusion
    // - where the extrusion ends

    isPatchEdge.setSize(mesh.nEdges());
    isPatchEdge = false;
    isPatchEdge.set(meshEdges);
    // Synchronise across coupled edges
    syncTools::syncEdgeList
    (
        mesh,
        isPatchEdge,
        orEqOp<unsigned int>(),
        false                   // initial value
    );

    isPatchBoundaryEdge.setSize(mesh.nEdges());
    isPatchBoundaryEdge = false;
    forAll(edgeGlobalFaces, edgei)
    {
        // Test that edge has single global extruded face.
        // Mark on processor that holds the face (since edgeGlobalFaces
        // only gets filled from pp faces so if there is only one this
        // is it)
        if (edgeGlobalFaces[edgei].size() == 1)
        {
            isPatchBoundaryEdge.set(meshEdges[edgei]);
        }
    }
    // Synchronise across coupled edges
    syncTools::syncEdgeList
    (
        mesh,
        isPatchBoundaryEdge,
        orEqOp<unsigned int>(),
        false                   // initial value
    );
}


void Foam::addPatchCellLayer::globalEdgeInfo
(
    const bool zoneFromAnyFace,

    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,

    labelList& patchEdgeToFace,     // face (in globalFaces index)
    labelList& patchEdgeToPatch,    // patch on face (or -1 for internal faces)
    labelList& patchEdgeToZone,     // zone on face
    bitSet& patchEdgeToFlip         // flip orientation on face
)
{
    // For every edge on the outside of the patch return a potential patch/
    // faceZone to extrude into.

    // Mark (mesh) edges on pp.
    bitSet isExtrudeEdge;
    bitSet isBoundaryEdge;
    markPatchEdges
    (
        mesh,
        pp,
        edgeGlobalFaces,
        meshEdges,

        isExtrudeEdge,
        isBoundaryEdge
    );

    // Build map of pp edges (in mesh point indexing). Note that this
    // is now also on processors that do not have pp (but do have the edge)
    EdgeMap<label> isBoundaryEdgeSet(pp.nEdges());
    for (const label edgei : isBoundaryEdge)
    {
        isBoundaryEdgeSet.insert(mesh.edges()[edgei], edgei);
    }
    EdgeMap<label> isExtrudeEdgeSet(pp.nEdges());
    for (const label edgei : isExtrudeEdge)
    {
        isExtrudeEdgeSet.insert(mesh.edges()[edgei], edgei);
    }


    bitSet isPpFace(mesh.nFaces());
    isPpFace.set(pp.addressing());
    // Note: no need to sync isPpFace since does not include processor patches


    const faceZoneMesh& fzs = mesh.faceZones();

    // Extract zone info into mesh face indexing for ease of addressing
    labelList faceToZone(mesh.nFaces(), -1);
    bitSet faceToFlip(mesh.nFaces());
    for (const faceZone& fz: fzs)
    {
        const labelList& addressing = fz;
        UIndirectList<label>(faceToZone, addressing) = fz.index();

        const boolList& fm = fz.flipMap();
        forAll(addressing, i)
        {
            faceToFlip[addressing[i]] = fm[i];
        }
    }


    // Storage (over all mesh edges)
    // - face that data originates from (in globalFaces indexing)
    labelList meshEdgeToFace(mesh.nEdges(), -1);
    // - patch (for boundary faces)
    labelList meshEdgeToPatch(mesh.nEdges(), -1);
    // - faceZone
    labelList meshEdgeToZone(mesh.nEdges(), -1);
    // - faceZone orientation
    bitSet meshEdgeToFlip(mesh.nEdges());

    //if (useInternalFaces)
    {
        const bitSet isInternalOrCoupled
        (
            syncTools::getInternalOrCoupledFaces(mesh)
        );

        // Loop over edges of the face to find any faceZone.
        // Edges kept as point pair so we don't construct mesh.faceEdges etc.

        for (const label facei : isInternalOrCoupled)
        {
            if (isPpFace[facei])
            {
                continue;
            }

            const face& f = mesh.faces()[facei];

            label prevPointi = f.last();
            for (const label pointi : f)
            {
                const edge e(prevPointi, pointi);

                // Check if edge is internal to extrusion. Take over faceZone
                // etc from internal face.
                const auto eFnd = isExtrudeEdgeSet.cfind(e);
                if (eFnd)
                {
                    const label edgei = eFnd();

                    if (faceToZone[facei] != -1)
                    {
                        // Found a zoned internal face. Use.
                        meshEdgeToFace[edgei] = globalFaces.toGlobal(facei);
                        meshEdgeToZone[edgei] = faceToZone[facei];
                        const edge& meshE = mesh.edges()[edgei];
                        const int d = edge::compare(e, meshE);
                        if (d == 1)
                        {
                            meshEdgeToFlip[edgei] = faceToFlip[facei];
                        }
                        else if (d == -1)
                        {
                            meshEdgeToFlip[edgei] = !faceToFlip[facei];
                        }
                        else
                        {
                            FatalErrorInFunction << "big problem"
                                << exit(FatalError);
                        }
                    }
                }

                prevPointi = pointi;
            }
        }
    }


    //if (useBoundaryFaces)
    {
        // Loop over all patches and find 'best' one (non-coupled,
        // non-extrusion, non-constraint(?)). Note that logic is slightly
        // different from internal faces loop above - first patch face
        // is being used instead of first zoned face for internal faces

        const polyBoundaryMesh& patches = mesh.boundaryMesh();


        for (const polyPatch& pp : patches)
        {
            if (!pp.coupled())
            {
                // TBD. Check for constraint? This is usually a geometric
                //      constraint and not a topological one so should
                //      be handled in the extrusion vector calculation instead?

                forAll(pp, i)
                {
                    const label facei = pp.start()+i;

                    if (!isPpFace[facei])
                    {
                        const face& f = pp[i];

                        label prevPointi = f.last();
                        for (const label pointi : f)
                        {
                            const edge e(prevPointi, pointi);
                            const auto eFnd =
                            (
                                zoneFromAnyFace
                              ? isExtrudeEdgeSet.cfind(e)
                              : isBoundaryEdgeSet.cfind(e)
                            );
                            if (eFnd)
                            {
                                const label edgei = eFnd();
                                if (meshEdgeToFace[edgei] == -1)
                                {
                                    // Found unassigned face. Use its
                                    // information.
                                    // Note that we use the lowest numbered
                                    // patch face.

                                    meshEdgeToFace[edgei] =
                                        globalFaces.toGlobal(facei);
                                }

                                // Override any patch info. Note that
                                // meshEdgeToFace might be an internal face.
                                if (meshEdgeToPatch[edgei] == -1)
                                {
                                    meshEdgeToPatch[edgei] = pp.index();
                                }

                                // Override any zone info
                                if (meshEdgeToZone[edgei] == -1)
                                {
                                    meshEdgeToZone[edgei] =
                                        faceToZone[facei];
                                    const edge& meshE = mesh.edges()[edgei];
                                    const int d = edge::compare(e, meshE);
                                    if (d == 1)
                                    {
                                        meshEdgeToFlip[edgei] =
                                            faceToFlip[facei];
                                    }
                                    else if (d == -1)
                                    {
                                        meshEdgeToFlip[edgei] =
                                           !faceToFlip[facei];
                                    }
                                    else
                                    {
                                        FatalErrorInFunction
                                            << "big problem"
                                            << exit(FatalError);
                                    }
                                }
                            }

                            prevPointi = pointi;
                        }
                    }
                }
            }
        }
    }


    // Synchronise across coupled edges. Max patch/face/faceZone wins
    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToFace,
        maxEqOp<label>(),
        label(-1)
    );
    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToPatch,
        maxEqOp<label>(),
        label(-1)
    );
    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToZone,
        maxEqOp<label>(),
        label(-1)
    );
//    // Note: flipMap not yet done. Needs edge orientation. This is handled
//    // later on.
//    if (Pstream::parRun())
//    {
//        const globalMeshData& gd = mesh.globalData();
//        const indirectPrimitivePatch& cpp = gd.coupledPatch();
//
//        labelList patchEdges;
//        labelList coupledEdges;
//        bitSet sameEdgeOrientation;
//        PatchTools::matchEdges
//        (
//            pp,
//            cpp,
//            patchEdges,
//            coupledEdges,
//            sameEdgeOrientation
//        );
//
//        // Convert data on pp edges to data on coupled patch
//        labelList cppEdgeZoneID(cpp.nEdges(), -1);
//        boolList cppEdgeFlip(cpp.nEdges(), false);
//        forAll(coupledEdges, i)
//        {
//            label cppEdgei = coupledEdges[i];
//            label ppEdgei = patchEdges[i];
//
//            cppEdgeZoneID[cppEdgei] = edgeZoneID[ppEdgei];
//            if (sameEdgeOrientation[i])
//            {
//                cppEdgeFlip[cppEdgei] = edgeFlip[ppEdgei];
//            }
//            else
//            {
//                cppEdgeFlip[cppEdgei] = !edgeFlip[ppEdgei];
//            }
//        }
//
//        // Sync
//        const globalIndexAndTransform& git = gd.globalTransforms();
//        const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();
//
//        globalMeshData::syncData
//        (
//            cppEdgeZoneID,
//            gd.globalEdgeSlaves(),
//            gd.globalEdgeTransformedSlaves(),
//            edgeMap,
//            git,
//            maxEqOp<label>(),
//            dummyTransform()
//        );
//        globalMeshData::syncData
//        (
//            cppEdgeFlip,
//            gd.globalEdgeSlaves(),
//            gd.globalEdgeTransformedSlaves(),
//            edgeMap,
//            git,
//            andEqOp<bool>(),
//            dummyTransform()
//        );
//
//        // Convert data on coupled edges to pp edges
//        forAll(coupledEdges, i)
//        {
//            label cppEdgei = coupledEdges[i];
//            label ppEdgei = patchEdges[i];
//
//            edgeZoneID[ppEdgei] = cppEdgeZoneID[cppEdgei];
//            if (sameEdgeOrientation[i])
//            {
//                edgeFlip[ppEdgei] = cppEdgeFlip[cppEdgei];
//            }
//            else
//            {
//                edgeFlip[ppEdgei] = !cppEdgeFlip[cppEdgei];
//            }
//        }
//    }

    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToFlip,
        orEqOp<unsigned int>(),
        0
    );

    // Extract pp info
    patchEdgeToFace = UIndirectList<label>(meshEdgeToFace, meshEdges);
    patchEdgeToPatch = UIndirectList<label>(meshEdgeToPatch, meshEdges);
    patchEdgeToZone = UIndirectList<label>(meshEdgeToZone, meshEdges);
    patchEdgeToFlip.setSize(meshEdges.size());
    patchEdgeToFlip = false;
    forAll(meshEdges, i)
    {
        patchEdgeToFlip[i] = meshEdgeToFlip[meshEdges[i]];
    }
}


void Foam::addPatchCellLayer::calcExtrudeInfo
(
    const bool zoneFromAnyFace,

    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const labelListList& globalEdgeFaces,
    const indirectPrimitivePatch& pp,

    labelList& edgePatchID,
    label& nPatches,
    Map<label>& nbrProcToPatch,
    Map<label>& patchToNbrProc,
    labelList& edgeZoneID,
    boolList& edgeFlip,
    labelList& inflateFaceID
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const globalMeshData& gd = mesh.globalData();

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    edgePatchID.setSize(pp.nEdges());
    edgePatchID = -1;
    nPatches = patches.size();
    nbrProcToPatch.clear();
    patchToNbrProc.clear();
    edgeZoneID.setSize(pp.nEdges());
    edgeZoneID = -1;
    edgeFlip.setSize(pp.nEdges());
    edgeFlip = false;
    inflateFaceID.setSize(pp.nEdges(), -1);


    // Determine properties for faces extruded from edges
    // - edge inbetween two different processors:
    //      - extrude as patch face on correct processor
    // - edge at end of patch (so edgeFaces.size() == 1):
    //      - use mesh face that is a boundary face
    // - edge internal to patch (so edgeFaces.size() == 2):


    // Pass1:
    // For all edges inbetween two processors: see if matches to existing
    // processor patch or create interprocessor-patch if necessary.
    // Sets edgePatchID[edgeI] but none of the other quantities

    forAll(globalEdgeFaces, edgei)
    {
        const labelList& eGlobalFaces = globalEdgeFaces[edgei];
        if
        (
            eGlobalFaces.size() == 2
         && pp.edgeFaces()[edgei].size() == 1
        )
        {
            // Locally but not globally a boundary edge. Hence a coupled
            // edge. Find the patch to use if on different processors.

            label f0 = eGlobalFaces[0];
            label f1 = eGlobalFaces[1];

            label otherProci = -1;
            if (globalFaces.isLocal(f0) && !globalFaces.isLocal(f1))
            {
                otherProci = globalFaces.whichProcID(f1);
            }
            else if (!globalFaces.isLocal(f0) && globalFaces.isLocal(f1))
            {
                otherProci = globalFaces.whichProcID(f0);
            }


            if (otherProci != -1)
            {
                // Use existing processorPolyPatch to otherProci?

                label procPatchi =
                    gd.topology().procPatchLookup(otherProci);

                if (procPatchi < 0)
                {
                    // No existing processorPolyPatch to otherProci.
                    // See if already marked for addition
                    procPatchi = nbrProcToPatch.lookup(otherProci, -1);

                    if (procPatchi < 0)
                    {
                        // Add new proc-patch, mark for addition.

                        procPatchi = nPatches;
                        ++nPatches;

                        nbrProcToPatch.insert(otherProci, procPatchi);
                        patchToNbrProc.insert(procPatchi, otherProci);
                    }
                }

                edgePatchID[edgei] = procPatchi;
            }
        }
    }


    // Pass2: determine face properties for all other edges
    // ----------------------------------------------------

    // Info for edges of pp
    labelList edgeToFace;
    labelList edgeToPatch;
    labelList edgeToZone;
    bitSet edgeToFlip;
    globalEdgeInfo
    (
        zoneFromAnyFace,    // internal edge info also from boundary faces

        mesh,
        globalFaces,
        globalEdgeFaces,
        pp,
        meshEdges,

        edgeToFace,         // face (in globalFaces index)
        edgeToPatch,        // patch on face (or -1 for internal faces)
        edgeToZone,         // zone on face
        edgeToFlip          // flip orientation on face
    );

    const labelListList& edgeFaces = pp.edgeFaces();

    DynamicList<label> dynMeshEdgeFaces;

    forAll(edgeFaces, edgei)
    {
        if (edgePatchID[edgei] == -1)
        {
            if (edgeFaces[edgei].size() == 2)
            {
                // Internal edge. Look at any face (internal or boundary) to
                // determine extrusion properties. First one that has zone
                // info wins
                if (globalFaces.isLocal(edgeToFace[edgei]))
                {
                    inflateFaceID[edgei] =
                        globalFaces.toLocal(edgeToFace[edgei]);
                }
                edgeZoneID[edgei] = edgeToZone[edgei];
                edgeFlip[edgei] = edgeToFlip[edgei];
            }
            else
            {
                // Proper, uncoupled patch edge. Boundary to get info from
                // might be on a different processor!

                if (globalFaces.isLocal(edgeToFace[edgei]))
                {
                    inflateFaceID[edgei] =
                        globalFaces.toLocal(edgeToFace[edgei]);
                }
                edgePatchID[edgei] = edgeToPatch[edgei];
                edgeZoneID[edgei] = edgeToZone[edgei];
                edgeFlip[edgei] = edgeToFlip[edgei];
            }
        }
    }



    // Now hopefully every boundary edge has a edge patch. Check
    if (debug)
    {
        forAll(edgeFaces, edgei)
        {
            if (edgeFaces[edgei].size() == 1 && edgePatchID[edgei] == -1)
            {
                const edge& e = pp.edges()[edgei];
                WarningInFunction
                    << "Have no sidePatchID for edge " << edgei << " points "
                    << pp.points()[pp.meshPoints()[e[0]]]
                    << pp.points()[pp.meshPoints()[e[1]]]
                    << endl;
            }
        }
    }



    // Pass3: for any faces set in pass1 see if we can find a processor face
    //        to inherit from (we only have a patch, not a patch face)
    forAll(edgeFaces, edgei)
    {
        if
        (
            edgeFaces[edgei].size() == 1
         && globalEdgeFaces[edgei].size() == 2
         && edgePatchID[edgei] != -1
         && inflateFaceID[edgei] == -1
        )
        {
            // 1. Do we have a local boundary face to inflate from

            label myFaceI = pp.addressing()[edgeFaces[edgei][0]];

            // Pick up any boundary face on this edge and use its properties
            label meshEdgei = meshEdges[edgei];
            const labelList& meshFaces = mesh.edgeFaces
            (
                meshEdgei,
                dynMeshEdgeFaces
            );

            forAll(meshFaces, k)
            {
                label facei = meshFaces[k];

                if (facei != myFaceI && !mesh.isInternalFace(facei))
                {
                    if (patches.whichPatch(facei) == edgePatchID[edgei])
                    {
                        setFaceProps
                        (
                            mesh,
                            pp,
                            edgei,
                            facei,

                            edgePatchID[edgei],
                            edgeZoneID[edgei],
                            edgeFlip[edgei],
                            inflateFaceID[edgei]
                        );
                        break;
                    }
                }
            }
        }
    }



    // Sync all data:
    // - edgePatchID : might be local in case of processor patch. Do not
    //   sync for now
    // - inflateFaceID: local. Do not sync
    // - edgeZoneID : global. sync.
    // - edgeFlip : global. sync.

    if (Pstream::parRun())
    {
        const globalMeshData& gd = mesh.globalData();
        const indirectPrimitivePatch& cpp = gd.coupledPatch();

        labelList patchEdges;
        labelList coupledEdges;
        bitSet sameEdgeOrientation;
        PatchTools::matchEdges
        (
            pp,
            cpp,
            patchEdges,
            coupledEdges,
            sameEdgeOrientation
        );

        // Convert data on pp edges to data on coupled patch
        labelList cppEdgeZoneID(cpp.nEdges(), -1);
        boolList cppEdgeFlip(cpp.nEdges(), false);
        forAll(coupledEdges, i)
        {
            label cppEdgei = coupledEdges[i];
            label ppEdgei = patchEdges[i];

            cppEdgeZoneID[cppEdgei] = edgeZoneID[ppEdgei];
            if (sameEdgeOrientation[i])
            {
                cppEdgeFlip[cppEdgei] = edgeFlip[ppEdgei];
            }
            else
            {
                cppEdgeFlip[cppEdgei] = !edgeFlip[ppEdgei];
            }
        }

        // Sync
        const globalIndexAndTransform& git = gd.globalTransforms();
        const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();

        globalMeshData::syncData
        (
            cppEdgeZoneID,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            edgeMap,
            git,
            maxEqOp<label>(),
            dummyTransform()
        );
        globalMeshData::syncData
        (
            cppEdgeFlip,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            edgeMap,
            git,
            andEqOp<bool>(),
            dummyTransform()
        );

        // Convert data on coupled edges to pp edges
        forAll(coupledEdges, i)
        {
            label cppEdgei = coupledEdges[i];
            label ppEdgei = patchEdges[i];

            edgeZoneID[ppEdgei] = cppEdgeZoneID[cppEdgei];
            if (sameEdgeOrientation[i])
            {
                edgeFlip[ppEdgei] = cppEdgeFlip[cppEdgei];
            }
            else
            {
                edgeFlip[ppEdgei] = !cppEdgeFlip[cppEdgei];
            }
        }
    }
}


void Foam::addPatchCellLayer::setRefinement
(
    const globalIndex& globalFaces,
    const labelListList& globalEdgeFaces,
    const scalarField& expansionRatio,
    const indirectPrimitivePatch& pp,
    const bitSet& ppFlip,

    const labelList& edgePatchID,
    const labelList& edgeZoneID,
    const boolList& edgeFlip,
    const labelList& inflateFaceID,

    const labelList& exposedPatchID,
    const labelList& nFaceLayers,
    const labelList& nPointLayers,
    const vectorField& firstLayerDisp,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "addPatchCellLayer::setRefinement : Adding up to "
            << gMax(nPointLayers)
            << " layers of cells to indirectPrimitivePatch with "
            << pp.nPoints() << " points" << endl;
    }

    if
    (
        pp.nPoints() != firstLayerDisp.size()
     || pp.nPoints() != nPointLayers.size()
     || pp.size() != nFaceLayers.size()
     || pp.size() != ppFlip.size()
    )
    {
        FatalErrorInFunction
            << "Size of new points is not same as number of points used by"
            << " the face subset" << endl
            << "  patch.nPoints:" << pp.nPoints()
            << "  displacement:" << firstLayerDisp.size()
            << "  nPointLayers:" << nPointLayers.size() << nl
            << " patch.nFaces:" << pp.size()
            << " flip map:" << ppFlip.size()
            << "  nFaceLayers:" << nFaceLayers.size()
            << abort(FatalError);
    }
    if (!addToMesh_)
    {
        // flip map should be false
        if (ppFlip.count())
        {
            FatalErrorInFunction
                << "In generating stand-alone mesh the flip map should be empty"
                << ". Instead it is " << ppFlip.count()
                << abort(FatalError);
        }
    }
    else if (debug)
    {
        // Maybe check for adding to neighbour of boundary faces? How about
        // coupled faces where the faceZone flipMap is negated

        // For all boundary faces:
        //  -1 : not extruded
        //   0 : extruded from owner outwards (flip = false)
        //   1 : extrude from neighbour outwards
        labelList stateAndFlip(mesh_.nBoundaryFaces(), 0);
        forAll(pp.addressing(), patchFacei)
        {
            if (nFaceLayers[patchFacei] > 0)
            {
                const label meshFacei = pp.addressing()[patchFacei];
                const label bFacei = meshFacei-mesh_.nInternalFaces();
                if (bFacei >= 0)
                {
                    stateAndFlip[bFacei] = label(ppFlip[patchFacei]);
                }
            }
        }
        // Make sure uncoupled patches do not trigger the error below
        for (const auto& patch : mesh_.boundaryMesh())
        {
            if (!patch.coupled())
            {
                forAll(patch, i)
                {
                    label& state = stateAndFlip[patch.offset()+i];
                    state = (state == 0 ? 1 : 0);
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, stateAndFlip);

        forAll(pp.addressing(), patchFacei)
        {
            if (nFaceLayers[patchFacei] > 0)
            {
                const label meshFacei = pp.addressing()[patchFacei];
                const label bFacei = meshFacei-mesh_.nInternalFaces();
                if (bFacei >= 0)
                {
                    if (stateAndFlip[bFacei] == -1)
                    {
                        FatalErrorInFunction
                            << "At extruded face:" << meshFacei
                            << " at:" << mesh_.faceCentres()[meshFacei]
                            << " locally have nLayers:"
                            << nFaceLayers[patchFacei]
                            << " but remotely have none" << exit(FatalError);
                    }
                    else if (stateAndFlip[bFacei] == label(ppFlip[patchFacei]))
                    {
                        FatalErrorInFunction
                            << "At extruded face:" << meshFacei
                            << " at:" << mesh_.faceCentres()[meshFacei]
                            << " locally have flip:" << ppFlip[patchFacei]
                            << " which is not the opposite of coupled version "
                            << stateAndFlip[bFacei]
                            << exit(FatalError);
                    }
                }
            }
        }
    }


    forAll(nPointLayers, i)
    {
        if (nPointLayers[i] < 0)
        {
            FatalErrorInFunction
                << "Illegal number of layers " << nPointLayers[i]
                << " at patch point " << i << abort(FatalError);
        }
    }
    forAll(nFaceLayers, i)
    {
        if (nFaceLayers[i] < 0)
        {
            FatalErrorInFunction
                << "Illegal number of layers " << nFaceLayers[i]
                << " at patch face " << i << abort(FatalError);
        }
    }

    forAll(globalEdgeFaces, edgei)
    {
        if (globalEdgeFaces[edgei].size() > 2)
        {
            const edge& e = pp.edges()[edgei];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                FatalErrorInFunction
                    << "Trying to extrude edge "
                    << e.line(pp.localPoints())
                    << " which is non-manifold (has "
                    << globalEdgeFaces[edgei].size()
                    << " local or coupled faces using it)"
                    << " of which "
                    << pp.edgeFaces()[edgei].size()
                    << " local"
                    << abort(FatalError);
            }
        }
    }


    const labelList& meshPoints = pp.meshPoints();


    // Determine which points are on which side of the extrusion
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    bitSet isBlockedFace(mesh_.nFaces());
    isBlockedFace.set(pp.addressing());

    // Some storage for edge-face-addressing.
    DynamicList<label> ef;

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh_.edges(), mesh_.pointEdges()));

    if (debug)
    {
        // Check synchronisation
        // ~~~~~~~~~~~~~~~~~~~~~

        {
            labelList n(mesh_.nPoints(), Zero);
            labelUIndList(n, meshPoints) = nPointLayers;
            syncTools::syncPointList(mesh_, n, maxEqOp<label>(), label(0));

            // Non-synced
            forAll(meshPoints, i)
            {
                label meshPointi = meshPoints[i];

                if (n[meshPointi] != nPointLayers[i])
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "On coupled point a different nLayers:"
                        << n[meshPointi] << " was specified."
                        << abort(FatalError);
                }
            }


            // Check that nPointLayers equals the max layers of connected faces
            // (or 0). Anything else makes no sense.
            labelList nFromFace(mesh_.nPoints(), Zero);
            forAll(nFaceLayers, i)
            {
                const face& f = pp[i];

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    nFromFace[pointi] = max(nFromFace[pointi], nFaceLayers[i]);
                }
            }
            syncTools::syncPointList
            (
                mesh_,
                nFromFace,
                maxEqOp<label>(),
                label(0)
            );

            forAll(nPointLayers, i)
            {
                label meshPointi = meshPoints[i];

                if
                (
                    nPointLayers[i] > 0
                 && nPointLayers[i] != nFromFace[meshPointi]
                )
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "but the max nLayers of surrounding faces is:"
                        << nFromFace[meshPointi]
                        << abort(FatalError);
                }
            }
        }

        {
            pointField d(mesh_.nPoints(), vector::max);
            UIndirectList<point>(d, meshPoints) = firstLayerDisp;
            syncTools::syncPointList
            (
                mesh_,
                d,
                minEqOp<vector>(),
                vector::max
            );

            forAll(meshPoints, i)
            {
                label meshPointi = meshPoints[i];

                if (mag(d[meshPointi] - firstLayerDisp[i]) > SMALL)
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified displacement:" << firstLayerDisp[i]
                        << endl
                        << "On coupled point a different displacement:"
                        << d[meshPointi] << " was specified."
                        << abort(FatalError);
                }
            }
        }

        // Check that edges of pp (so ones that become boundary faces)
        // connect to only one boundary face. Guarantees uniqueness of
        // patch that they go into so if this is a coupled patch both
        // sides decide the same.
        // ~~~~~~~~~~~~~~~~~~~~~~

        for (label edgei = pp.nInternalEdges(); edgei < pp.nEdges(); edgei++)
        {
            const edge& e = pp.edges()[edgei];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                // Edge is to become a face

                const labelList& eFaces = pp.edgeFaces()[edgei];

                // First check: pp should be single connected.
                if (eFaces.size() != 1)
                {
                    FatalErrorInFunction
                        << "boundary-edge-to-be-extruded:"
                        << pp.points()[meshPoints[e[0]]]
                        << pp.points()[meshPoints[e[1]]]
                        << " has more than two faces using it:" << eFaces
                        << abort(FatalError);
                }

                //label myFacei = pp.addressing()[eFaces[0]];
                //
                //label meshEdgei = meshEdges[edgei];
                //
                //// Mesh faces using edge
                //const labelList& meshFaces = mesh_.edgeFaces(meshEdgei, ef);
                //
                //// Check that there is only one patchface using edge.
                //const polyBoundaryMesh& patches = mesh_.boundaryMesh();
                //
                //label bFacei = -1;
                //
                //forAll(meshFaces, i)
                //{
                //    label facei = meshFaces[i];
                //
                //    if (facei != myFacei)
                //    {
                //        if (!mesh_.isInternalFace(facei))
                //        {
                //            if (bFacei == -1)
                //            {
                //                bFacei = facei;
                //            }
                //            else
                //            {
                //                //FatalErrorInFunction
                //                WarningInFunction
                //                    << "boundary-edge-to-be-extruded:"
                //                    << pp.points()[meshPoints[e[0]]]
                //                    << pp.points()[meshPoints[e[1]]]
                //                    << " has more than one boundary faces"
                //                    << " using it:"
                //                    << bFacei << " fc:"
                //                    << mesh_.faceCentres()[bFacei]
                //                    << " patch:" << patches.whichPatch(bFacei)
                //                    << " and " << facei << " fc:"
                //                    << mesh_.faceCentres()[facei]
                //                    << " patch:" << patches.whichPatch(facei)
                //                    << endl;
                //                    //abort(FatalError);
                //            }
                //        }
                //    }
                //}
            }
        }
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Precalculated patchID for each patch face
    labelList patchID(pp.size());

    forAll(pp, patchFacei)
    {
        label meshFacei = pp.addressing()[patchFacei];

        patchID[patchFacei] = patches.whichPatch(meshFacei);
    }


    // From master point (in patch point label) to added points (in mesh point
    // label)
    addedPoints_.setSize(pp.nPoints());

    // Mark points that do not get extruded by setting size of addedPoints_ to 0
    label nTruncated = 0;

    forAll(nPointLayers, patchPointi)
    {
        if (nPointLayers[patchPointi] > 0)
        {
            addedPoints_[patchPointi].setSize(nPointLayers[patchPointi]);
        }
        else
        {
            nTruncated++;
        }
    }

    if (debug)
    {
        Pout<< "Not adding points at " << nTruncated << " out of "
            << pp.nPoints() << " points" << endl;
    }


    // Store per face whether it uses the duplicated point or the original one
    // Also mark any affected cells. We could transport the duplicated point
    // itself but since it is a processor-local index only we only transport
    // a boolean.
    // Per face, per index in face either -1 or a valid index. Note:
    // most faces are not affected in which case the face will be zero size
    // and only have a nullptr and a size. Note: we do not actually use the
    // value of the index - it might be a remote patch point so not make any
    // sense. All we want to know is whether it is -1 or not.
    faceList baseFaces(mesh_.nFaces());
    {
        bitSet isDupPatchPoint(pp.nPoints());
        forAll(nPointLayers, patchPointi)
        {
            if (nPointLayers[patchPointi] > 0)
            {
                isDupPatchPoint.set(patchPointi);
            }
        }
        // Knock out point if not all local faces using it are extruded
        // Guess not needed since done outside of routine already.
        forAll(nFaceLayers, patchFacei)
        {
            if (nFaceLayers[patchFacei] == 0)
            {
                const face& f = pp.localFaces()[patchFacei];
                for (const label patchPointi : f)
                {
                    if (isDupPatchPoint[patchPointi])
                    {
                        isDupPatchPoint.unset(patchPointi);
                    }
                }
            }
        }

        findDuplicatedPoints
        (
            mesh_,
            pp,
            ppFlip,           // optional orientation on top of pp
            isBlockedFace,    // any mesh faces not to be traversed.
                              // Usually pp.addressing()
            isDupPatchPoint,
            extrude_,

            baseFaces
        );
    }


    //
    // Create new points
    //

    // If creating new mesh: copy existing patch points
    labelList copiedPatchPoints;
    if (!addToMesh_)
    {
        copiedPatchPoints.setSize(firstLayerDisp.size());
        forAll(firstLayerDisp, patchPointi)
        {
            if (addedPoints_[patchPointi].size())
            {
                label meshPointi = meshPoints[patchPointi];
                label zoneI = mesh_.pointZones().whichZone(meshPointi);
                copiedPatchPoints[patchPointi] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        mesh_.points()[meshPointi],         // point
                        -1,         // master point
                        zoneI,      // zone for point
                        true        // supports a cell
                    )
                );
            }
        }
    }


    // Create points for additional layers
    forAll(firstLayerDisp, patchPointi)
    {
        if (addedPoints_[patchPointi].size())
        {
            const label meshPointi = meshPoints[patchPointi];
            const label zoneI = mesh_.pointZones().whichZone(meshPointi);

            point pt = mesh_.points()[meshPointi];

            vector disp = firstLayerDisp[patchPointi];

            forAll(addedPoints_[patchPointi], i)
            {
                pt += disp;

                const label addedVertI = meshMod.setAction
                (
                    polyAddPoint
                    (
                        pt,         // point
                        (addToMesh_ ? meshPointi : -1), // master point
                        zoneI,      // zone for point
                        true        // supports a cell
                    )
                );


                //Pout<< "Adding point:" << addedVertI << " at:" << pt
                //    << " from point:" << meshPointi
                //    << " at:" << mesh_.points()[meshPointi]
                //    << endl;

                addedPoints_[patchPointi][i] = addedVertI;

                disp *= expansionRatio[patchPointi];
            }
        }
    }


    //
    // Add cells to all boundaryFaces
    //

    labelListList addedCells(pp.size());

    forAll(pp, patchFacei)
    {
        if (nFaceLayers[patchFacei] > 0)
        {
            const label meshFacei = pp.addressing()[patchFacei];

            label extrudeCelli = -2;
            label extrudeZonei;
            if (!addToMesh_)
            {
                extrudeCelli = -1;
                const label ownCelli = mesh_.faceOwner()[meshFacei];
                extrudeZonei = mesh_.cellZones().whichZone(ownCelli);
            }
            else if (!ppFlip[patchFacei])
            {
                // Normal: extrude from owner side
                extrudeCelli = mesh_.faceOwner()[meshFacei];
                extrudeZonei = mesh_.cellZones().whichZone(extrudeCelli);
            }
            else if (mesh_.isInternalFace(meshFacei))
            {
                // Extrude from neighbour side (if internal). Might be
                // that it is a coupled face and the other side is
                // extruded
                extrudeCelli = mesh_.faceNeighbour()[meshFacei];
                extrudeZonei = mesh_.cellZones().whichZone(extrudeCelli);
            }

            if (extrudeCelli != -2)
            {
                addedCells[patchFacei].setSize(nFaceLayers[patchFacei]);

                for (label i = 0; i < nFaceLayers[patchFacei]; i++)
                {
                    // Note: add from cell (owner of patch face) or from face?
                    // for now add from cell so we can map easily.
                    addedCells[patchFacei][i] = meshMod.setAction
                    (
                        polyAddCell
                        (
                            -1,             // master point
                            -1,             // master edge
                            -1,             // master face
                            extrudeCelli,   // master cell
                            extrudeZonei    // zone for cell
                        )
                    );

                    //Pout<< "Added cell:" << addedCells[patchFacei][i]
                    //    << " from master:" << extrudeCelli
                    //    << " at:" << mesh_.cellCentres()[extrudeCelli]
                    //    << endl;
                }
            }
        }
    }



    // Create faces from the original patch faces.
    // These faces are created from original patch faces outwards (or inwards)
    // so in order
    // of increasing cell number. So orientation should be same as original
    // patch face for them to have owner<neighbour.

    layerFaces_.setSize(pp.size());

    forAll(pp.localFaces(), patchFacei)
    {
        label meshFacei = pp.addressing()[patchFacei];

        if (addedCells[patchFacei].size())
        {
            layerFaces_[patchFacei].setSize(addedCells[patchFacei].size() + 1);

            // Get duplicated vertices on the patch face.
            const face& f = pp.localFaces()[patchFacei];

            face newFace(f.size());

            forAll(addedCells[patchFacei], i)
            {
                forAll(f, fp)
                {
                    if (addedPoints_[f[fp]].empty())
                    {
                        // Keep original point
                        newFace[fp] =
                        (
                            addToMesh_
                          ? meshPoints[f[fp]]
                          : copiedPatchPoints[f[fp]]
                        );
                    }
                    else
                    {
                        // Get new outside point
                        label offset =
                            addedPoints_[f[fp]].size()
                          - addedCells[patchFacei].size();
                        newFace[fp] = addedPoints_[f[fp]][i+offset];
                    }
                }
                //Pout<< "   newFace:" << newFace << endl;
                //Pout<< "   coords:"
                //    << UIndirectList<point>(meshMod.points(), newFace)
                //    << " normal:" << newFace.unitNormal(meshMod.points())
                //    << endl;

                // Get new neighbour
                label own;
                label nei = -1;
                label patchi = -1;
                label zonei = -1;
                bool zoneFlip = false;
                bool fluxFlip = false;

                if (i == addedCells[patchFacei].size()-1)
                {
                    // Last layer so
                    //  - extrude : original patch or other cell
                    //  - intrude : original cell

                    if (extrude_)
                    {
                        patchi = patchID[patchFacei];
                        zonei = mesh_.faceZones().whichZone(meshFacei);
                        if (zonei != -1)
                        {
                            const faceZone& fz = mesh_.faceZones()[zonei];
                            zoneFlip = fz.flipMap()[fz.whichFace(meshFacei)];
                        }
                        if (patchi == -1)
                        {
                            // Internal face between added cell and original
                            // neighbour / owner
                            own =
                            (
                               !ppFlip[patchFacei]
                              ? mesh_.faceNeighbour()[meshFacei]
                              : mesh_.faceOwner()[meshFacei]
                            );
                            nei = addedCells[patchFacei][i];
                            if (!ppFlip[patchFacei])
                            {
                                newFace = newFace.reverseFace();
                                if (zonei != -1)
                                {
                                    zoneFlip = !zoneFlip;
                                }
                                fluxFlip = true;
                            }
                        }
                        else
                        {
                            // Boundary face
                            own = addedCells[patchFacei][i];
                        }
                    }
                    else
                    {
                        // Intrude : face connects to original cell
                        //if (patchID[patchFacei] == -1)
                        {
                            // Internal face between added cell and original
                            // owner / neighbour. TBD: what if ppFlip set on
                            // boundary face (so with no neighbour)
                            own =
                            (
                               !ppFlip[patchFacei]
                              ? mesh_.faceOwner()[meshFacei]
                              : mesh_.faceNeighbour()[meshFacei]
                            );
                            nei = addedCells[patchFacei][i];
                            if (ppFlip[patchFacei])
                            {
                                // Since neighbour now owns the face flip it
                                newFace = newFace.reverseFace();
                                if (zonei != -1)
                                {
                                    zoneFlip = !zoneFlip;
                                }
                                fluxFlip = true;
                            }
                        }
                    }
                }
                else
                {
                    // Internal face between layer i and i+1
                    own = addedCells[patchFacei][i];
                    nei = addedCells[patchFacei][i+1];

                    if (ppFlip[patchFacei])
                    {
                        // Since neighbour now owns the face flip it
                        newFace = newFace.reverseFace();
                        fluxFlip = true;
                    }
                }

                layerFaces_[patchFacei][i+1] = meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,                    // face
                        own,                        // owner
                        nei,                        // neighbour
                        -1,                         // master point
                        -1,                         // master edge
                        (addToMesh_ ? meshFacei : -1), // master face
                        fluxFlip,                   // flux flip
                        patchi,                     // patch for face
                        zonei,                      // zone for face
                        zoneFlip                    // face zone flip
                    )
                );

                //Pout<< "added layer face:" << layerFaces_[patchFacei][i+1]
                //    << " verts:" << newFace
                //    << " newFc:" << newFace.centre(meshMod.points())
                //    << " originalFc:" << mesh_.faceCentres()[meshFacei]
                //    << nl
                //    << "     n:" << newFace.unitNormal(meshMod.points())
                //    << " own:" << own << " nei:" << nei
                //    << " patchi:" << patchi
                //    << " zonei:" << zonei
                //    << endl;
            }
        }
    }

    //
    // Modify original face to have addedCells as neighbour
    //

    if (addToMesh_)
    {
        forAll(pp, patchFacei)
        {
            if (addedCells[patchFacei].size())
            {
                const label meshFacei = pp.addressing()[patchFacei];

                layerFaces_[patchFacei][0] = meshFacei;
                const face& f = pp[patchFacei];

                if (extrude_)
                {
                    const label own =
                    (
                        !ppFlip[patchFacei]
                      ? mesh_.faceOwner()[meshFacei]
                      : mesh_.faceNeighbour()[meshFacei]
                    );
                    const label nei = addedCells[patchFacei][0];

                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            (ppFlip[patchFacei] ? f.reverseFace() : f),// verts
                            meshFacei,                      // label of face
                            own,                            // owner
                            nei,                            // neighbour
                            ppFlip[patchFacei],             // face flip
                            -1,                             // patch for face
                            true, //false,                  // remove from zone
                            -1, //zoneI,                    // zone for face
                            false                           // face flip in zone
                        )
                    );
                }
                else
                {
                    label patchi;
                    label zonei;
                    bool zoneOrient;
                    setFaceProps
                    (
                        mesh_,
                        meshFacei,

                        patchi,
                        zonei,
                        zoneOrient
                    );

                    label own;
                    label nei = -1;
                    bool flip = false;
                    if (!ppFlip[patchFacei])
                    {
                        if (patchi == -1)
                        {
                            own = mesh_.faceNeighbour()[meshFacei];
                            nei = addedCells[patchFacei].first();
                            flip = true;
                        }
                        else
                        {
                            own = addedCells[patchFacei].first();
                        }
                    }
                    else
                    {
                        if (patchi == -1)
                        {
                            own = mesh_.faceOwner()[meshFacei];
                            nei = addedCells[patchFacei].first();
                        }
                        else
                        {
                            own = addedCells[patchFacei].first();
                        }
                    }

                    if (zonei != -1 && flip)
                    {
                        zoneOrient = !zoneOrient;
                    }

                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            (flip ? f.reverseFace() : f),   // verts
                            meshFacei,                      // label of face
                            own,                            // owner
                            nei,                            // neighbour
                            flip,                           // face flip
                            patchi,                         // patch for face
                            false,                          // remove from zone
                            zonei,                          // zone for face
                            zoneOrient                      // face flip in zone
                        )
                    );
                }

                //Pout<< "Modified original face " << meshFacei
                //    << " at:" << mesh_.faceCentres()[meshFacei]
                //    << " new own:" << meshMod.faceOwner()[meshFacei]
                //    << " new nei:" << meshMod.faceNeighbour()[meshFacei]
                //    << " verts:" << meshMod.faces()[meshFacei]
                //    << " n:"
                //    << meshMod.faces()[meshFacei].unitNormal(meshMod.points())
                //    << endl;
            }
        }
    }
    else
    {
        // If creating new mesh: reverse original faces and put them
        // in the exposed patch ID.
        forAll(pp, patchFacei)
        {
            if (addedCells[patchFacei].size())
            {
                label meshFacei = pp.addressing()[patchFacei];
                label zoneI = mesh_.faceZones().whichZone(meshFacei);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    const faceZone& fz = mesh_.faceZones()[zoneI];
                    zoneFlip = !fz.flipMap()[fz.whichFace(meshFacei)];
                }

                // Reverse and renumber old patch face.
                face f(pp.localFaces()[patchFacei].reverseFace());
                forAll(f, fp)
                {
                    f[fp] = copiedPatchPoints[f[fp]];
                }

                layerFaces_[patchFacei][0] = meshMod.setAction
                (
                    polyAddFace
                    (
                        f,                          // modified face
                        addedCells[patchFacei][0],  // owner
                        -1,                         // neighbour
                        -1,                         // masterPoint
                        -1,                         // masterEdge
                        -1,                         // masterFace
                        true,                       // face flip
                        exposedPatchID[patchFacei], // patch for face
                        zoneI,                      // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
        }
    }



    //
    // Create 'side' faces, one per edge that is being extended.
    //

    const labelListList& faceEdges = pp.faceEdges();
    const faceList& localFaces = pp.localFaces();
    const edgeList& edges = pp.edges();

    // Get number of layers per edge. This is 0 if edge is not extruded;
    // max of connected faces otherwise.
    labelList edgeLayers(pp.nEdges());

    {
        // Use list over mesh.nEdges() since syncTools does not yet support
        // partial list synchronisation.
        labelList meshEdgeLayers(mesh_.nEdges(), -1);

        forAll(meshEdges, edgei)
        {
            const edge& e = edges[edgei];

            label meshEdgei = meshEdges[edgei];

            if ((nPointLayers[e[0]] == 0) && (nPointLayers[e[1]] == 0))
            {
                meshEdgeLayers[meshEdgei] = 0;
            }
            else
            {
                const labelList& eFaces = pp.edgeFaces()[edgei];

                forAll(eFaces, i)
                {
                    meshEdgeLayers[meshEdgei] = max
                    (
                        nFaceLayers[eFaces[i]],
                        meshEdgeLayers[meshEdgei]
                    );
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            meshEdgeLayers,
            maxEqOp<label>(),
            label(0)            // initial value
        );

        forAll(meshEdges, edgei)
        {
            edgeLayers[edgei] = meshEdgeLayers[meshEdges[edgei]];
        }
    }


    // Mark off which edges have been extruded
    boolList doneEdge(pp.nEdges(), false);


    // Create faces. Per face walk connected edges and find string of edges
    // between the same two faces and extrude string into a single face.
    forAll(pp, patchFacei)
    {
        // Get edges of face in vertex order
        const labelList& fEdges = faceEdges[patchFacei];

        forAll(fEdges, i)
        {
            // Get string of edges that needs to be extruded as a single face.
            // Returned as indices in fEdges.
            labelPair indexPair
            (
                getEdgeString
                (
                    pp,
                    globalEdgeFaces,
                    doneEdge,
                    patchFacei,
                    globalFaces.toGlobal(pp.addressing()[patchFacei])
                )
            );

            //Pout<< "Found unextruded edges in edges:" << fEdges
            //    << " start:" << indexPair[0]
            //    << " end:" << indexPair[1]
            //    << endl;

            const label startFp = indexPair[0];
            const label endFp = indexPair[1];

            if (startFp != -1)
            {
                // Extrude edges from indexPair[0] up to indexPair[1]
                // (note indexPair = indices of edges. There is one more vertex
                //  than edges)
                const face& f = localFaces[patchFacei];

                labelList stringedVerts;
                if (endFp >= startFp)
                {
                    stringedVerts.setSize(endFp-startFp+2);
                }
                else
                {
                    stringedVerts.setSize(endFp+f.size()-startFp+2);
                }

                label fp = startFp;

                for (label i = 0; i < stringedVerts.size()-1; i++)
                {
                    stringedVerts[i] = f[fp];
                    doneEdge[fEdges[fp]] = true;
                    fp = f.fcIndex(fp);
                }
                stringedVerts.last() = f[fp];


                // Now stringedVerts contains the vertices in order of face f.
                // This is consistent with the order if f becomes the owner cell
                // and nbrFacei the neighbour cell. Note that the cells get
                // added in order of pp so we can just use face ordering and
                // because we loop in incrementing order as well we will
                // always have nbrFacei > patchFacei.

                const label startEdgei = fEdges[startFp];

                const label meshEdgei = meshEdges[startEdgei];

                const label numEdgeSideFaces = edgeLayers[startEdgei];

                for (label i = 0; i < numEdgeSideFaces; i++)
                {
                    const label vEnd = stringedVerts.last();
                    const label vStart = stringedVerts[0];

                    // calculate number of points making up a face
                    label newFp = 2*stringedVerts.size();

                    if (i == 0)
                    {
                        // layer 0 gets all the truncation of neighbouring
                        // faces with more layers.
                        if (addedPoints_[vEnd].size())
                        {
                            newFp +=
                                addedPoints_[vEnd].size() - numEdgeSideFaces;
                        }
                        if (addedPoints_[vStart].size())
                        {
                            newFp +=
                                addedPoints_[vStart].size() - numEdgeSideFaces;
                        }
                    }

                    face newFace(newFp);

                    newFp = 0;

                    // For layer 0 get pp points, for all other layers get
                    // points of layer-1.
                    if (i == 0)
                    {
                        forAll(stringedVerts, stringedI)
                        {
                            const label v = stringedVerts[stringedI];
                            addVertex
                            (
                                (
                                    addToMesh_
                                  ? meshPoints[v]
                                  : copiedPatchPoints[v]
                                ),
                                newFace,
                                newFp
                            );
                        }
                    }
                    else
                    {
                        forAll(stringedVerts, stringedI)
                        {
                            const label v = stringedVerts[stringedI];
                            if (addedPoints_[v].size())
                            {
                                const label offset =
                                    addedPoints_[v].size() - numEdgeSideFaces;
                                addVertex
                                (
                                    addedPoints_[v][i+offset-1],
                                    newFace,
                                    newFp
                                );
                            }
                            else
                            {
                                addVertex
                                (
                                    (
                                        addToMesh_
                                      ? meshPoints[v]
                                      : copiedPatchPoints[v]
                                    ),
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    // add points between stringed vertices (end)
                    if (numEdgeSideFaces < addedPoints_[vEnd].size())
                    {
                        if (i == 0 && addedPoints_[vEnd].size())
                        {
                            const label offset =
                                addedPoints_[vEnd].size() - numEdgeSideFaces;
                            for (label ioff = 0; ioff < offset; ioff++)
                            {
                                addVertex
                                (
                                    addedPoints_[vEnd][ioff],
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    forAllReverse(stringedVerts, stringedI)
                    {
                        const label v = stringedVerts[stringedI];
                        if (addedPoints_[v].size())
                        {
                            const label offset =
                                addedPoints_[v].size() - numEdgeSideFaces;
                            addVertex
                            (
                                addedPoints_[v][i+offset],
                                newFace,
                                newFp
                            );
                        }
                        else
                        {
                            addVertex
                            (
                                (
                                    addToMesh_
                                  ? meshPoints[v]
                                  : copiedPatchPoints[v]
                                ),
                                newFace,
                                newFp
                            );
                        }
                    }


                    // add points between stringed vertices (start)
                    if (numEdgeSideFaces < addedPoints_[vStart].size())
                    {
                        if (i == 0 && addedPoints_[vStart].size())
                        {
                            const label offset =
                                addedPoints_[vStart].size() - numEdgeSideFaces;
                            for (label ioff = offset-1; ioff >= 0; ioff--)
                            {
                                addVertex
                                (
                                    addedPoints_[vStart][ioff],
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    if (newFp >= 3)
                    {
                        // Add face inbetween faces patchFacei and nbrFacei
                        // (possibly -1 for external edges)

                        newFace.setSize(newFp);

                        // Walked edges as if owner face was extruded. Reverse
                        // for neighbour face extrusion.
                        if (extrude_ == ppFlip[patchFacei])
                        {
                            newFace = newFace.reverseFace();
                        }

                        if (debug)
                        {
                            labelHashSet verts(2*newFace.size());
                            forAll(newFace, fp)
                            {
                                if (!verts.insert(newFace[fp]))
                                {
                                    FatalErrorInFunction
                                        << "Duplicate vertex in face"
                                        << " to be added." << nl
                                        << "newFace:" << newFace << nl
                                        << "points:"
                                        <<  UIndirectList<point>
                                            (
                                                meshMod.points(),
                                                newFace
                                            ) << nl
                                        << "Layer:" << i
                                        << " out of:" << numEdgeSideFaces << nl
                                        << "ExtrudeEdge:" << meshEdgei
                                        << " at:"
                                        <<  mesh_.edges()[meshEdgei].line
                                            (
                                                mesh_.points()
                                            ) << nl
                                        << "string:" << stringedVerts
                                        << "stringpoints:"
                                        << UIndirectList<point>
                                            (
                                                pp.localPoints(),
                                                stringedVerts
                                            ) << nl
                                        << "stringNLayers:"
                                        <<  labelUIndList
                                            (
                                                nPointLayers,
                                                stringedVerts
                                            ) << nl
                                        << abort(FatalError);
                                }
                            }
                        }

                        label nbrFacei = nbrFace
                        (
                            pp.edgeFaces(),
                            startEdgei,
                            patchFacei
                        );

                        const labelList& meshFaces = mesh_.edgeFaces
                        (
                            meshEdgei,
                            ef
                        );

                        // Because we walk in order of patch face and in order
                        // of face edges so face orientation will be opposite
                        // that of the patch edge
                        bool zoneFlip = false;
                        if (edgeZoneID[startEdgei] != -1)
                        {
                            zoneFlip = !edgeFlip[startEdgei];
                        }

                        addSideFace
                        (
                            pp,
                            addedCells,

                            newFace,                // vertices of new face
                            edgePatchID[startEdgei],// -1 or patch for face
                            edgeZoneID[startEdgei],
                            zoneFlip,
                            inflateFaceID[startEdgei],

                            patchFacei,
                            nbrFacei,
                            meshEdgei,          // (mesh) edge to inflate
                            i,                  // layer
                            numEdgeSideFaces,   // num layers
                            meshFaces,          // edgeFaces
                            meshMod
                        );
                    }
                }
            }
        }
    }


    // Adjust side faces if they're on the side where points were duplicated
    // (i.e. adding to internal faces). Like duplicatePoints::setRefinement.
    if (addToMesh_)
    {
        face newFace;

        // For intrusion correction later on: store displacement of modified
        // points so they can be adapted on connected faces.
        pointField displacement(mesh_.nPoints(), Zero);

        forAll(mesh_.faces(), facei)
        {
            const face& f = mesh_.faces()[facei];
            const face& baseF = baseFaces[facei];

            if (isBlockedFace(facei) || baseF.empty())
            {
                // Either part of patch or no duplicated points on face
                continue;
            }

            // Start off from original face
            newFace = f;
            forAll(f, fp)
            {
                const label meshPointi = f[fp];
                if (baseF[fp] != -1)
                {
                    // Duplicated point. Might not be present on processor?
                    //const label patchPointi = pp.meshPointMap()[meshPointi];
                    const auto pointFnd = pp.meshPointMap().find(meshPointi);
                    if (pointFnd)
                    {
                        const label patchPointi = pointFnd();
                        const label addedPointi =
                            addedPoints_[patchPointi].last();

                        // Adapt vertices of face
                        newFace[fp] = addedPointi;

                        // Store displacement for syncing later on
                        displacement[meshPointi] =
                            meshMod.points()[addedPointi]
                           -mesh_.points()[meshPointi];
                    }
                }
            }

            //Pout<< "for face:" << facei << nl
            //    << "    old:" << f
            //    << " n:" << newFace.unitNormal(meshMod.points())
            //    << " coords:" << UIndirectList<point>(mesh_.points(), f)
            //    << nl
            //    << "    new:" << newFace
            //    << " n:" << newFace.unitNormal(meshMod.points())
            //    << " coords:"
            //    << UIndirectList<point>(meshMod.points(), newFace)
            //    << endl;

            // Get current zone info
            label zoneID = mesh_.faceZones().whichZone(facei);
            bool zoneFlip = false;
            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }


            if (mesh_.isInternalFace(facei))
            {
                meshMod.modifyFace
                (
                    newFace,                    // modified face
                    facei,                      // label of modified face
                    mesh_.faceOwner()[facei],   // owner
                    mesh_.faceNeighbour()[facei],   // neighbour
                    false,                      // face flip
                    -1,                         // patch for face
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                );
            }
            else
            {
                meshMod.modifyFace
                (
                    newFace,                    // modified face
                    facei,                      // label of modified face
                    mesh_.faceOwner()[facei],   // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    patches.whichPatch(facei),  // patch for face
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                );
            }


            //Pout<< "Adapted point on existing face:" << facei
            //    << " at:" << mesh_.faceCentres()[facei]
            //    << " zone:" << zoneID << nl
            //    << "    old:" << f
            //    << " n:" << newFace.unitNormal(meshMod.points())
            //    << " coords:" << UIndirectList<point>(mesh_.points(), f)
            //    << nl
            //    << "    new:" << newFace
            //    << " n:" << newFace.unitNormal(meshMod.points())
            //    << " coords:"
            //    << UIndirectList<point>(meshMod.points(), newFace)
            //    << endl;
        }

        if (!extrude_)
        {
            // Bit tricky:
            // - if intruding we're modifying the vertices on existing
            //   coupled faces (in extrude mode we're adding additional points/
            //   faces where we're already doing the right thing)
            // - if the other side face is itself not on an extrudePatch
            //   it will not be adapted so you'll get a coupled point mismatch
            // - so send over the new point position (or use some map on
            //   coupled faces to describe which vertex on the face is changed?)
            // - and adapt any point that was marked but not an extrusion point
            // - to test: take 2x2 blockMesh and remove one of the cells. Put
            //   all 3 cells on different processors. Now one of the processors
            //   will be coupled to extrude/intruded points but itself have
            //   no patch.

            syncTools::syncPointList
            (
                mesh_,
                displacement,
                maxMagSqrEqOp<vector>(),
                vector::zero
            );

            forAll(baseFaces, facei)
            {
                const face& f = mesh_.faces()[facei];
                const face& baseF = baseFaces[facei];

                if (isBlockedFace(facei) || baseF.empty())
                {
                    // Either part of patch or no duplicated points on face
                    continue;
                }

                // Start off from original face
                forAll(f, fp)
                {
                    if (baseF[fp] != -1)
                    {
                        // Duplicated point. Might not be present on processor?
                        const label meshPointi = f[fp];

                        if
                        (
                           !pp.meshPointMap().found(meshPointi)
                         && displacement[meshPointi] != vector::zero
                        )
                        {
                            const point newPt
                            (
                                mesh_.points()[meshPointi]
                              + displacement[meshPointi]
                            );
                            meshMod.modifyPoint
                            (
                                meshPointi,
                                newPt,
                                mesh_.pointZones().whichZone(meshPointi),
                                true
                            );

                            // Unmark as being done
                            displacement[meshPointi] = Zero;
                        }
                    }
                }
            }
        }
    }



    //if (debug & 4)
    //{
    //    Pout<< "Checking whole mesh" << endl;
    //    pointField fCtrs;
    //    vectorField fAreas;
    //    const UList<face>& faces = meshMod.faces();
    //    const pointField p(meshMod.points());
    //    primitiveMeshTools::makeFaceCentresAndAreas
    //    (
    //        faces,
    //        p,
    //        fCtrs,
    //        fAreas
    //    );
    //
    //    const label nCells = max(meshMod.faceOwner())+1;
    //
    //    pointField cellCtrs;
    //    scalarField cellVols;
    //    {
    //        // See primitiveMeshTools::makeCellCentresAndVols
    //        // Clear the fields for accumulation
    //        cellCtrs.setSize(nCells, Zero);
    //        cellVols.setSize(nCells, Zero);
    //
    //        const labelList& own = meshMod.faceOwner();
    //        const labelList& nei = meshMod.faceNeighbour();
    //
    //        Field<solveVector> cEst(nCells, Zero);
    //        labelField nCellFaces(nCells, Zero);
    //
    //        forAll(own, facei)
    //        {
    //            cEst[own[facei]] += solveVector(fCtrs[facei]);
    //            ++nCellFaces[own[facei]];
    //        }
    //
    //        forAll(nei, facei)
    //        {
    //            if (nei[facei] == -1)
    //            {
    //                continue;
    //            }
    //            cEst[nei[facei]] += solveVector(fCtrs[facei]);
    //            ++nCellFaces[nei[facei]];
    //        }
    //
    //        forAll(cEst, celli)
    //        {
    //            cEst[celli] /= nCellFaces[celli];
    //        }
    //
    //        forAll(own, facei)
    //        {
    //            const solveVector fc(fCtrs[facei]);
    //            const solveVector fA(fAreas[facei]);
    //
    //            // Calculate 3*face-pyramid volume
    //            solveScalar pyr3Vol = fA & (fc - cEst[own[facei]]);
    //
    //            // Calculate face-pyramid centre
    //            solveVector pc = (3.0/4.0)*fc + (1.0/4.0)*cEst[own[facei]];
    //
    //            // Accumulate volume-weighted face-pyramid centre
    //            cellCtrs[own[facei]] += pyr3Vol*pc;
    //
    //            // Accumulate face-pyramid volume
    //            cellVols[own[facei]] += pyr3Vol;
    //        }
    //
    //        forAll(nei, facei)
    //        {
    //            if (nei[facei] == -1)
    //            {
    //                continue;
    //            }
    //
    //            const solveVector fc(fCtrs[facei]);
    //            const solveVector fA(fAreas[facei]);
    //
    //            // Calculate 3*face-pyramid volume
    //            solveScalar pyr3Vol = fA & (cEst[nei[facei]] - fc);
    //
    //            // Calculate face-pyramid centre
    //            solveVector pc = (3.0/4.0)*fc + (1.0/4.0)*cEst[nei[facei]];
    //
    //            // Accumulate volume-weighted face-pyramid centre
    //            cellCtrs[nei[facei]] += pyr3Vol*pc;
    //
    //            // Accumulate face-pyramid volume
    //            cellVols[nei[facei]] += pyr3Vol;
    //        }
    //
    //        forAll(cellCtrs, celli)
    //        {
    //            if (mag(cellVols[celli]) > VSMALL)
    //            {
    //                cellCtrs[celli] /= cellVols[celli];
    //            }
    //            else
    //            {
    //                cellCtrs[celli] = cEst[celli];
    //            }
    //        }
    //        cellVols *= (1.0/3.0);
    //    }
    //    cellList cells(nCells);
    //    {
    //        const labelList& own = meshMod.faceOwner();
    //        const labelList& nei = meshMod.faceNeighbour();
    //
    //        labelList nFaces(nCells, 0);
    //        forAll(own, facei)
    //        {
    //            nFaces[own[facei]]++;
    //        }
    //        forAll(nei, facei)
    //        {
    //            if (nei[facei] == -1)
    //            {
    //                continue;
    //            }
    //            nFaces[nei[facei]]++;
    //        }
    //
    //        forAll(cells, celli)
    //        {
    //            cells[celli].resize_nocopy(nFaces[celli]);
    //            nFaces[celli] = 0;
    //        }
    //
    //        forAll(own, facei)
    //        {
    //            const label celli = own[facei];
    //            cells[celli][nFaces[celli]++] = facei;
    //        }
    //        forAll(nei, facei)
    //        {
    //            if (nei[facei] == -1)
    //            {
    //                continue;
    //            }
    //            const label celli = nei[facei];
    //            cells[celli][nFaces[celli]++] = facei;
    //        }
    //    }
    //
    //
    //    scalarField openness;
    //    {
    //        const labelList& own = meshMod.faceOwner();
    //        const labelList& nei = meshMod.faceNeighbour();
    //
    //        // Loop through cell faces and sum up the face area vectors for
    //        // each cell. This should be zero in all vector components
    //
    //        vectorField sumClosed(nCells, Zero);
    //        vectorField sumMagClosed(nCells, Zero);
    //
    //        forAll(own, facei)
    //        {
    //            // Add to owner
    //            sumClosed[own[facei]] += fAreas[facei];
    //            sumMagClosed[own[facei]] += cmptMag(fAreas[facei]);
    //        }
    //
    //        forAll(nei, facei)
    //        {
    //            // Subtract from neighbour
    //            if (nei[facei] == -1)
    //            {
    //                continue;
    //            }
    //
    //            sumClosed[nei[facei]] -= fAreas[facei];
    //            sumMagClosed[nei[facei]] += cmptMag(fAreas[facei]);
    //        }
    //
    //
    //        // Check the sums
    //        openness.setSize(nCells);
    //
    //        forAll(sumClosed, celli)
    //        {
    //            scalar maxOpenness = 0;
    //
    //            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    //            {
    //                maxOpenness = max
    //                (
    //                    maxOpenness,
    //                    mag(sumClosed[celli][cmpt])
    //                   /(sumMagClosed[celli][cmpt] + ROOTVSMALL)
    //                );
    //            }
    //            openness[celli] = maxOpenness;
    //
    //            if (openness[celli] > 1e-9)
    //            {
    //                Pout<< "cell:" << celli
    //                    << " at:" << cellCtrs[celli]
    //                    << " openness:" << openness[celli]
    //                    << endl;
    //                OBJstream os
    //                (
    //                    mesh_.time().timePath()
    //                   /"cell_" + Foam::name(celli) + ".obj"
    //                );
    //                faceList cellFaces
    //                (
    //                    UIndirectList<face>
    //                    (
    //                        meshMod.faces(),
    //                        cells[celli]
    //                    )
    //                );
    //                os.write(cellFaces, mesh_.points(), false);
    //                Pout<< "Written " << os.nVertices() << " to "
    //                    << os.name() << endl;
    //            }
    //        }
    //    }
    //}
}


void Foam::addPatchCellLayer::updateMesh
(
    const mapPolyMesh& morphMap,
    const labelList& faceMap,   // new to old patch faces
    const labelList& pointMap   // new to old patch points
)
{
    {
        labelListList newAddedPoints(pointMap.size());

        forAll(newAddedPoints, newPointi)
        {
            label oldPointi = pointMap[newPointi];

            const labelList& added = addedPoints_[oldPointi];

            labelList& newAdded = newAddedPoints[newPointi];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newPointi = morphMap.reversePointMap()[added[i]];

                if (newPointi >= 0)
                {
                    newAdded[newI++] = newPointi;
                }
            }
            newAdded.setSize(newI);
        }
        addedPoints_.transfer(newAddedPoints);
    }

    {
        labelListList newLayerFaces(faceMap.size());

        forAll(newLayerFaces, newFacei)
        {
            label oldFacei = faceMap[newFacei];

            const labelList& added = layerFaces_[oldFacei];

            labelList& newAdded = newLayerFaces[newFacei];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newFacei = morphMap.reverseFaceMap()[added[i]];

                if (newFacei >= 0)
                {
                    newAdded[newI++] = newFacei;
                }
            }
            newAdded.setSize(newI);
        }
        layerFaces_.transfer(newLayerFaces);
    }
}


void Foam::addPatchCellLayer::findDuplicatedPoints
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const bitSet& ppFlip,           // optional orientation on top of pp
    const bitSet& isBlockedFace,    // any mesh faces not to be traversed.
                                    // Usually pp.addressing()
    const bitSet& isDupPatchPoint,
    const bool extrude,
    faceList& isDupMeshPoint    // per face, per index either -1
                                // or some index (value not relevant)
)
{
    if (isBlockedFace.size() != mesh.nFaces())
    {
        FatalErrorInFunction << "Incorrect size" << exit(FatalError);
    }
    if (ppFlip.size() != pp.size() || isDupPatchPoint.size() != pp.nPoints())
    {
        FatalErrorInFunction << "Incorrect patch sizes"
            << exit(FatalError);
    }

    const auto& owner = mesh.faceOwner();
    const auto& neighbour = mesh.faceNeighbour();


    // Store per face whether it uses the duplicated point or the original one
    // Also mark any affected cells. We could transport the duplicated point
    // itself but since it is a processor-local index only we only transport
    // a boolean.

    // Per face, per index in face either -1 or a valid index. Note:
    // most faces are not affected in which case the face will be zero size
    // and only have a nullptr and a size.
    isDupMeshPoint.resize_nocopy(mesh.nFaces());
    isDupMeshPoint = face();
    bitSet isAffectedCell(mesh.nCells());
    {
        const faceList& localFaces = pp.localFaces();
        forAll(localFaces, patchFacei)
        {
            const face& f = localFaces[patchFacei];
            forAll(f, fp)
            {
                const label patchPointi = f[fp];
                if (isDupPatchPoint[patchPointi])
                {
                    const label meshFacei = pp.addressing()[patchFacei];
                    face& baseF = isDupMeshPoint[meshFacei];
                    // Initialise to -1 if not yet sized
                    baseF.setSize(f.size(), -1);
                    baseF[fp] = pp.meshPoints()[patchPointi];

                    if (extrude == ppFlip[patchFacei])
                    {
                        // either:
                        // - extrude out of neighbour so new points connect
                        //   to owner
                        // - or intrude into owner
                        isAffectedCell.set(owner[meshFacei]);
                    }
                    else if (mesh.isInternalFace(meshFacei))
                    {
                        // Owner unaffected. Affected points on neighbour side
                        isAffectedCell.set(neighbour[meshFacei]);
                    }
                }
            }
        }
    }


    // Transport affected side across faces. Could do across edges: say we have
    // a loose cell edge-(but not face-)connected to face-to-be-extruded do
    // we want it to move with the extrusion or stay connected to the original?
    // For now just keep it connected to the original.
    {
        // Work space
        Map<label> minPointValue;

        while (true)
        {
            bitSet newIsAffectedCell(mesh.nCells());

            label nChanged = 0;
            for (const label celli : isAffectedCell)
            {
                const cell& cFaces = mesh.cells()[celli];

                // 1. Determine marked base points. Inside a single cell all
                //    faces (e.g. 3 for hex) use the same 'instance' of a point.
                minPointValue.clear();
                for (const label facei : cFaces)
                {
                    const face& baseF = isDupMeshPoint[facei];
                    const face& f = mesh.faces()[facei];

                    forAll(baseF, fp)
                    {
                        if (baseF[fp] != -1)
                        {
                            // Could check here for inconsistent patchPoint
                            // e.g. cell using both sides of a
                            // face-to-be-extruded.

                            const auto mpm = pp.meshPointMap().find(f[fp]);
                            if (mpm && !isDupPatchPoint[mpm()])
                            {
                                // Local copy of point is explicitly not
                                // marked for extrusion so ignore. This
                                // occasionally happens with point-connected
                                // faces where one face (& point) does
                                // extrude but the other face does not
                                // since pointing in the other direction.
                                // Should ideally be covered by checking
                                // across edges, not across points.
                            }
                            else
                            {
                                minPointValue.insert(f[fp], baseF[fp]);
                            }
                        }
                    }
                }


                // 2. Transport marked points on all cell points
                for (const label facei : cFaces)
                {
                    const face& f = mesh.faces()[facei];
                    face& baseF = isDupMeshPoint[facei];

                    const label oldNChanged = nChanged;
                    forAll(f, fp)
                    {
                        const auto fnd = minPointValue.find(f[fp]);
                        if (fnd.good())
                        {
                            const auto mpm = pp.meshPointMap().find(f[fp]);
                            if (mpm && !isDupPatchPoint[mpm()])
                            {
                                // See above
                                continue;
                            }

                            baseF.setSize(f.size(), -1);
                            if (baseF[fp] == -1)
                            {
                                baseF[fp] = fnd();
                                nChanged++;
                            }
                        }
                    }

                    if (!isBlockedFace(facei) && nChanged > oldNChanged)
                    {
                        // Mark neighbouring cells. Note that we do NOT check
                        // for whether cell is already in isAffectedCell but
                        // always add it. This is because different information
                        // can come in from the neighbour in different
                        // iterations.
                        newIsAffectedCell.set(owner[facei]);
                        if (mesh.isInternalFace(facei))
                        {
                            newIsAffectedCell.set(neighbour[facei]);
                        }
                    }
                }
            }


            if (debug)
            {
                Pout<< "isAffectedCell:" << isAffectedCell.count() << endl;
                Pout<< "newIsAffectedCell:" << newIsAffectedCell.count()
                    << endl;
                Pout<< "nChanged:" << nChanged << endl;
            }


            if (!returnReduceOr(nChanged))
            {
                break;
            }


            // Transport minimum across coupled faces
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Faces can be either
            // - null
            // - face with some -1 or >= 0
            // - blocked face (isBlockedFace, should be synced already)
            faceList l(mesh.nBoundaryFaces());
            forAll(l, bFacei)
            {
                const label facei = mesh.nInternalFaces()+bFacei;

                if (isBlockedFace(facei))
                {
                    // Make sure nothing is transferred. Since isBlockedFace is
                    // synced both sides should have null
                    if (l[bFacei].size())
                    {
                        l[bFacei].clear();
                    }
                }
                else
                {
                    l[bFacei] = isDupMeshPoint[facei];
                }
            }


            syncTools::syncBoundaryFaceList
            (
                mesh,
                l,
                combineEqOpFace(),      // transport vertex >= 0
                Foam::dummyTransform()  // dummy transformation
            );

            // Copy back
            forAll(l, bFacei)
            {
                const label facei = mesh.nInternalFaces()+bFacei;

                if (!isBlockedFace(facei))
                {
                    // 1. Check if anything changed. Update isAffectedCell.
                    // Note: avoid special handling of comparing zero-sized
                    //       faces (see face::operator==). Review.
                    const labelUList& newVts = l[bFacei];
                    const labelUList& oldVts = isDupMeshPoint[facei];

                    if (newVts != oldVts)
                    {
                        const label own = owner[facei];
                        if (!isAffectedCell[own])
                        {
                            newIsAffectedCell.set(own);
                        }
                    }


                    // 2. Update isDupMeshPoint
                    isDupMeshPoint[facei] = l[bFacei];
                }
            }


            isAffectedCell = newIsAffectedCell;
        }
    }
}


// ************************************************************************* //
