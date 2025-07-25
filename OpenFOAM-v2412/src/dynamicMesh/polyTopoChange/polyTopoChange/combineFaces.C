/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "combineFaces.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyRemoveFace.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "polyAddPoint.H"
#include "syncTools.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combineFaces, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Test face for (almost) convexeness. Allows certain convexity before
// complaining.
// For every two consecutive edges calculate the normal. If it differs too
// much from the average face normal complain.
bool Foam::combineFaces::convexFace
(
    const scalar minConcaveCos,
    const pointField& points,
    const face& f
)
{
    // Get outwards pointing normal of f, only the sign matters.
    const vector areaNorm = f.areaNormal(points);

    // Normalized vector from f[size-1] to f[0];
    vector ePrev(points[f.first()] - points[f.last()]);
    scalar magEPrev = mag(ePrev);
    ePrev /= magEPrev + VSMALL;

    forAll(f, fp0)
    {
        // Normalized vector between two consecutive points
        vector e10(points[f.nextLabel(fp0)] - points[f.thisLabel(fp0)]);
        scalar magE10 = mag(e10);
        e10 /= magE10 + VSMALL;

        if (magEPrev > SMALL && magE10 > SMALL)
        {
            vector edgeNormal = ePrev ^ e10;

            if ((edgeNormal & areaNorm) < 0)
            {
                // Concave. Check angle.
                if ((ePrev & e10) < minConcaveCos)
                {
                    return false;
                }
            }
        }

        ePrev = e10;
        magEPrev = magE10;
    }

    // Not a single internal angle is concave so face is convex.
    return true;
}


// Determines if set of faces is valid to collapse into single face.
bool Foam::combineFaces::validFace
(
    const scalar minConcaveCos,
    const indirectPrimitivePatch& bigFace
)
{
    // Get outside vertices (in local vertex numbering)
    const labelListList& edgeLoops = bigFace.edgeLoops();

    if (edgeLoops.size() > 1)
    {
        return false;
    }

    bool isNonManifold = bigFace.checkPointManifold(false, nullptr);
    if (isNonManifold)
    {
        return false;
    }

    // Check for convexness
    face f(getOutsideFace(bigFace));

    return convexFace(minConcaveCos, bigFace.points(), f);
}


void Foam::combineFaces::regioniseFaces
(
    const scalar minCos,
    const bool mergeAcrossPatches,
    const label celli,
    const labelList& cEdges,
    Map<label>& faceRegion
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (const label edgei : cEdges)
    {
        label f0, f1;
        meshTools::getEdgeFaces(mesh_, celli, edgei, f0, f1);

        const vector& a0 = mesh_.faceAreas()[f0];
        const vector& a1 = mesh_.faceAreas()[f1];

        const label p0 = patches.whichPatch(f0);
        const label p1 = patches.whichPatch(f1);

        // Face can be merged if
        // - small angle
        // - mergeAcrossPatches=false : same non-constraint patch
        // - mergeAcrossPatches=true  : always (if non-constraint patch)
        //   (this logic could be extended to e.g. merge faces on symm plane
        //    if they have similar normals. But there might be lots of other
        //    constraints which disallow merging so this decision ideally should
        //    be up to patch type)
        if
        (
            p0 != -1
         && p1 != -1
         && !(
                polyPatch::constraintType(patches[p0].type())
             || polyPatch::constraintType(patches[p1].type())
            )
        )
        {
            if (!mergeAcrossPatches && (p0 != p1))
            {
                continue;
            }

            const vector f0Normal = normalised(a0);
            const vector f1Normal = normalised(a1);

            if ((f0Normal & f1Normal) > minCos)
            {
                const label region0 = faceRegion.lookup(f0, -1);
                const label region1 = faceRegion.lookup(f1, -1);

                if (region0 == -1)
                {
                    if (region1 == -1)
                    {
                        const label useRegion = faceRegion.size();
                        faceRegion.insert(f0, useRegion);
                        faceRegion.insert(f1, useRegion);
                    }
                    else
                    {
                        faceRegion.insert(f0, region1);
                    }
                }
                else
                {
                    if (region1 == -1)
                    {
                        faceRegion.insert(f1, region0);
                    }
                    else if (region0 != region1)
                    {
                        // Merge the two regions
                        const label useRegion = min(region0, region1);
                        const label freeRegion = max(region0, region1);

                        forAllIters(faceRegion, iter)
                        {
                            if (iter.val() == freeRegion)
                            {
                                iter.val() = useRegion;
                            }
                        }
                    }
                }
            }
        }
    }
}


bool Foam::combineFaces::faceNeighboursValid
(
    const label celli,
    const Map<label>& faceRegion
) const
{
    if (faceRegion.size() <= 1)
    {
        return true;
    }

    const cell& cFaces = mesh_.cells()[celli];

    DynamicList<label> storage;

    // Test for face collapsing to edge since too many neighbours merged.
    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        if (!faceRegion.found(facei))
        {
            const labelList& fEdges = mesh_.faceEdges(facei, storage);

            // Count number of remaining faces neighbouring facei. This has
            // to be 3 or more.

            // Unregioned neighbouring faces
            DynamicList<label> neighbourFaces(cFaces.size());
            // Regioned neighbouring faces
            labelHashSet neighbourRegions;

            forAll(fEdges, i)
            {
                const label edgeI = fEdges[i];
                label nbrI = meshTools::otherFace(mesh_, celli, facei, edgeI);

                const auto iter = faceRegion.cfind(nbrI);

                if (iter.good())
                {
                    neighbourRegions.insert(iter.val());
                }
                else
                {
                    neighbourFaces.push_uniq(nbrI);
                }
            }

            if ((neighbourFaces.size()+neighbourRegions.size()) < 3)
            {
                return false;
            }
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::combineFaces::combineFaces
(
    const polyMesh& mesh,
    const bool undoable
)
:
    mesh_(mesh),
    undoable_(undoable),
    masterFace_(0),
    faceSetsVertices_(0),
    savedPointLabels_(0),
    savedPoints_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::combineFaces::getMergeSets
(
    const scalar featureCos,
    const scalar minConcaveCos,
    const labelHashSet& boundaryCells,
    const bool mergeAcrossPatches
) const
{
    // Lists of faces that can be merged.
    DynamicList<labelList> allFaceSets(boundaryCells.size() / 10);
    // Storage for on-the-fly cell-edge addressing.
    labelHashSet set;
    DynamicList<label> storage;

    // On all cells regionise the faces
    for (const label celli : boundaryCells)
    {
        const cell& cFaces = mesh_.cells()[celli];

        const labelList& cEdges = mesh_.cellEdges(celli, set, storage);

        // Region per face
        Map<label> faceRegion(cFaces.size());
        regioniseFaces
        (
            featureCos,
            mergeAcrossPatches,
            celli,
            cEdges,
            faceRegion
        );

        // Now we have in faceRegion for every face the region with planar
        // face sharing the same region. We now check whether the resulting
        // sets cause a face
        // - to become a set of edges since too many faces are merged.
        // - to become convex

        if (faceNeighboursValid(celli, faceRegion))
        {
            // Create region-to-faces addressing
            Map<labelList> regionToFaces(faceRegion.size());

            forAllConstIters(faceRegion, iter)
            {
                const label facei = iter.key();
                const label region = iter.val();

                auto regionFnd = regionToFaces.find(region);

                if (regionFnd.good())
                {
                    labelList& setFaces = regionFnd();
                    label sz = setFaces.size();
                    setFaces.setSize(sz+1);
                    setFaces[sz] = facei;
                }
                else
                {
                    regionToFaces.insert(region, labelList(1, facei));
                }
            }

            // For every set check if it forms a valid convex face

            forAllIters(regionToFaces, iter)
            {
                // Make face out of setFaces
                indirectPrimitivePatch bigFace
                (
                    IndirectList<face>
                    (
                        mesh_.faces(),
                        iter.val()
                    ),
                    mesh_.points()
                );

                // Only store if -only one outside loop -which forms convex face
                if (validFace(minConcaveCos, bigFace))
                {
                    labelList& faceIDs = iter.val();

                    // For cross-patch merging we want to make the
                    // largest face the one to decide the final patch
                    // (i.e. master face)
                    if (mergeAcrossPatches)
                    {
                        const vectorField& areas = mesh_.faceAreas();

                        label maxIndex = 0;
                        scalar maxMagSqr = magSqr(areas[faceIDs[0]]);
                        for (label i = 1; i < faceIDs.size(); ++i)
                        {
                            const scalar a2 = magSqr(areas[faceIDs[i]]);
                            if (a2 > maxMagSqr)
                            {
                                maxMagSqr = a2;
                                maxIndex = i;
                            }
                        }
                        if (maxIndex != 0)
                        {
                            std::swap(faceIDs[0], faceIDs[maxIndex]);
                        }
                    }

                    allFaceSets.append(faceIDs);
                }
            }
        }
    }

    return allFaceSets.shrink();
}


Foam::labelListList Foam::combineFaces::getMergeSets
(
    const scalar featureCos,
    const scalar minConcaveCos,
    const bool mergeAcrossPatches
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Pick up all cells on boundary
    labelHashSet boundaryCells(mesh_.nBoundaryFaces());

    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];

        if (!patch.coupled())
        {
            forAll(patch, i)
            {
                boundaryCells.insert(mesh_.faceOwner()[patch.start()+i]);
            }
        }
    }

    return getMergeSets
    (
        featureCos,
        minConcaveCos,
        boundaryCells,
        mergeAcrossPatches
    );
}


// Gets outside edgeloop as a face
// - in same order as faces
// - in mesh vertex labels
Foam::face Foam::combineFaces::getOutsideFace
(
    const indirectPrimitivePatch& fp
)
{
    if (fp.edgeLoops().size() != 1)
    {
        FatalErrorInFunction
            << "Multiple outside loops:" << fp.edgeLoops()
            << abort(FatalError);
    }

    // Get first boundary edge. Since guaranteed one edgeLoop when in here this
    // edge must be on it.
    label bEdgeI = fp.nInternalEdges();

    const edge& e = fp.edges()[bEdgeI];

    const labelList& eFaces = fp.edgeFaces()[bEdgeI];

    if (eFaces.size() != 1)
    {
        FatalErrorInFunction
            << "boundary edge:" << bEdgeI
            << " points:" << fp.meshPoints()[e[0]]
            << ' ' << fp.meshPoints()[e[1]]
            << " on indirectPrimitivePatch has " << eFaces.size()
            << " faces using it" << abort(FatalError);
    }


    // Outside loop
    const labelList& outsideLoop = fp.edgeLoops()[0];


    // Get order of edge e in outside loop
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // edgeLoopConsistent : true  = edge in same order as outsideloop
    //                      false = opposite order
    bool edgeLoopConsistent = false;

    {
        label index0 = outsideLoop.find(e[0]);
        label index1 = outsideLoop.find(e[1]);

        if (index0 == -1 || index1 == -1)
        {
            FatalErrorInFunction
                << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in edgeLoop:" << outsideLoop << abort(FatalError);
        }
        else if (index1 == outsideLoop.fcIndex(index0))
        {
            edgeLoopConsistent = true;
        }
        else if (index0 == outsideLoop.fcIndex(index1))
        {
            edgeLoopConsistent = false;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " on consecutive points in edgeLoop:"
                << outsideLoop << abort(FatalError);
        }
    }


    // Get face in local vertices
    const face& localF = fp.localFaces()[eFaces[0]];

    // Get order of edge in localF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // faceEdgeConsistent : true  = edge in same order as localF
    //                      false = opposite order
    bool faceEdgeConsistent = false;

    {
        // Find edge in face.
        label index = fp.faceEdges()[eFaces[0]].find(bEdgeI);

        if (index == -1)
        {
            FatalErrorInFunction
                << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in face:" << eFaces[0]
                << " edges:" << fp.faceEdges()[eFaces[0]]
                << abort(FatalError);
        }
        else if (localF[index] == e[0] && localF.nextLabel(index) == e[1])
        {
            faceEdgeConsistent = true;
        }
        else if (localF[index] == e[1] && localF.nextLabel(index) == e[0])
        {
            faceEdgeConsistent = false;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in face:" << eFaces[0] << " verts:" << localF
                << abort(FatalError);
        }
    }

    // Get face in mesh points.
    face meshFace(renumber(fp.meshPoints(), outsideLoop));

    if (faceEdgeConsistent != edgeLoopConsistent)
    {
        reverse(meshFace);
    }
    return meshFace;
}


void Foam::combineFaces::setRefinement
(
    const labelListList& faceSets,
    polyTopoChange& meshMod
)
{
    if (undoable_)
    {
        masterFace_.setSize(faceSets.size());
        faceSetsVertices_.setSize(faceSets.size());
        forAll(faceSets, setI)
        {
            const labelList& setFaces = faceSets[setI];

            masterFace_[setI] = setFaces[0];
            faceSetsVertices_[setI] = UIndirectList<face>
            (
                mesh_.faces(),
                setFaces
            );
        }
    }

    // Running count of number of faces using a point
    labelList nPointFaces(mesh_.nPoints(), Zero);

    const labelListList& pointFaces = mesh_.pointFaces();

    forAll(pointFaces, pointi)
    {
        nPointFaces[pointi] = pointFaces[pointi].size();
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(faceSets, setI)
    {
        const labelList& setFaces = faceSets[setI];

        // Check
        if (debug)
        {
            forAll(setFaces, i)
            {
                label patchi = patches.whichPatch(setFaces[i]);

                if (patchi == -1 || patches[patchi].coupled())
                {
                    FatalErrorInFunction
                        << "Can only merge non-coupled boundary faces"
                        << " but found internal or coupled face:"
                        << setFaces[i] << " in set " << setI
                        << abort(FatalError);
                }
            }
        }

        // Make face out of setFaces
        indirectPrimitivePatch bigFace
        (
            IndirectList<face>
            (
                mesh_.faces(),
                setFaces
            ),
            mesh_.points()
        );

        // Get outside vertices (in local vertex numbering)
        const labelListList& edgeLoops = bigFace.edgeLoops();

        if (edgeLoops.size() != 1)
        {
            FatalErrorInFunction
                << "Faces to-be-merged " << setFaces
                << " do not form a single big face." << nl
                << abort(FatalError);
        }


        // Delete all faces except master, whose face gets modified.

        // Modify master face
        // ~~~~~~~~~~~~~~~~~~

        label masterFacei = setFaces[0];

        // Get outside face in mesh vertex labels
        face outsideFace(getOutsideFace(bigFace));

        label zoneID = mesh_.faceZones().whichZone(masterFacei);

        bool zoneFlip = false;

        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];

            zoneFlip = fZone.flipMap()[fZone.whichFace(masterFacei)];
        }

        label patchi = mesh_.boundaryMesh().whichPatch(masterFacei);

        meshMod.setAction
        (
            polyModifyFace
            (
                outsideFace,                    // modified face
                masterFacei,                    // label of face being modified
                mesh_.faceOwner()[masterFacei], // owner
                -1,                             // neighbour
                false,                          // face flip
                patchi,                         // patch for face
                false,                          // remove from zone
                zoneID,                         // zone for face
                zoneFlip                        // face flip in zone
            )
        );


        // Delete all non-master faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label i = 1; i < setFaces.size(); i++)
        {
            meshMod.setAction(polyRemoveFace(setFaces[i]));
        }


        // Mark unused points for deletion
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Remove count of all faces combined
        forAll(setFaces, i)
        {
            const face& f = mesh_.faces()[setFaces[i]];

            forAll(f, fp)
            {
                nPointFaces[f[fp]]--;
            }
        }
        // But keep count on new outside face
        forAll(outsideFace, fp)
        {
            nPointFaces[outsideFace[fp]]++;
        }
    }


    // If one or more processors use the point it needs to be kept.
    syncTools::syncPointList
    (
        mesh_,
        nPointFaces,
        plusEqOp<label>(),
        label(0)            // null value
    );

    // Remove all unused points. Store position if undoable.
    if (!undoable_)
    {
        forAll(nPointFaces, pointi)
        {
            if (nPointFaces[pointi] == 0)
            {
                meshMod.setAction(polyRemovePoint(pointi));
            }
        }
    }
    else
    {
        // Count removed points
        label n = 0;
        forAll(nPointFaces, pointi)
        {
            if (nPointFaces[pointi] == 0)
            {
                n++;
            }
        }

        savedPointLabels_.setSize(n);
        savedPoints_.setSize(n);
        Map<label> meshToSaved(2*n);

        // Remove points and store position
        n = 0;
        forAll(nPointFaces, pointi)
        {
            if (nPointFaces[pointi] == 0)
            {
                meshMod.setAction(polyRemovePoint(pointi));

                savedPointLabels_[n] = pointi;
                savedPoints_[n] = mesh_.points()[pointi];

                meshToSaved.insert(pointi, n);
                n++;
            }
        }

        // Update stored vertex labels. Negative indices index into local points
        forAll(faceSetsVertices_, setI)
        {
            faceList& setFaces = faceSetsVertices_[setI];

            forAll(setFaces, i)
            {
                face& f = setFaces[i];

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    if (nPointFaces[pointi] == 0)
                    {
                        f[fp] = -meshToSaved[pointi]-1;
                    }
                }
            }
        }
    }
}


void Foam::combineFaces::updateMesh(const mapPolyMesh& map)
{
    if (undoable_)
    {
        // Master face just renumbering of point labels
        inplaceRenumber(map.reverseFaceMap(), masterFace_);

        // Stored faces refer to backed-up vertices (not changed)
        // and normal vertices (need to be renumbered)
        forAll(faceSetsVertices_, setI)
        {
            faceList& faces = faceSetsVertices_[setI];

            forAll(faces, i)
            {
                // Note: faces[i] can have negative elements.
                face& f = faces[i];

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    if (pointi >= 0)
                    {
                        f[fp] = map.reversePointMap()[pointi];

                        if (f[fp] < 0)
                        {
                            FatalErrorInFunction
                                << "In set " << setI << " at position " << i
                                << " with master face "
                                << masterFace_[setI] << nl
                                << "the points of the slave face " << faces[i]
                                << " don't exist anymore."
                                << abort(FatalError);
                        }
                    }
                }
            }
        }
    }
}


// Note that faces on coupled patches are never combined (or at least
// getMergeSets never picks them up. Thus the points on them will never get
// deleted since still used by the coupled faces. So never a risk of undoing
// things (faces/points) on coupled patches.
void Foam::combineFaces::setUnrefinement
(
    const labelList& masterFaces,
    polyTopoChange& meshMod,
    Map<label>& restoredPoints,
    Map<label>& restoredFaces,
    Map<label>& restoredCells
)
{
    if (!undoable_)
    {
        FatalErrorInFunction
            << "Can only call setUnrefinement if constructed with"
            << " unrefinement capability." << exit(FatalError);
    }


    // Restored points
    labelList addedPoints(savedPoints_.size(), -1);

    // Invert set-to-master-face
    Map<label> masterToSet(masterFace_.size());

    forAll(masterFace_, setI)
    {
        if (masterFace_[setI] >= 0)
        {
            masterToSet.insert(masterFace_[setI], setI);
        }
    }

    forAll(masterFaces, i)
    {
        const label masterFacei = masterFaces[i];

        const auto iter = masterToSet.cfind(masterFacei);

        if (!iter.good())
        {
            FatalErrorInFunction
                << "Master face " << masterFacei
                << " is not the master of one of the merge sets"
                << " or has already been merged"
                << abort(FatalError);
        }

        const label setI = iter.val();


        // Update faces of the merge setI for reintroduced vertices
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        faceList& faces = faceSetsVertices_[setI];

        if (faces.empty())
        {
            FatalErrorInFunction
                << "Set " << setI << " with master face " << masterFacei
                << " has already been merged." << abort(FatalError);
        }

        forAll(faces, i)
        {
            face& f = faces[i];

            forAll(f, fp)
            {
                label pointi = f[fp];

                if (pointi < 0)
                {
                    label localI = -pointi-1;

                    if (addedPoints[localI] == -1)
                    {
                        // First occurrence of saved point. Reintroduce point
                        addedPoints[localI] = meshMod.setAction
                        (
                            polyAddPoint
                            (
                                savedPoints_[localI],   // point
                                -1,                     // master point
                                -1,                     // zone for point
                                true                    // supports a cell
                            )
                        );
                        restoredPoints.insert
                        (
                            addedPoints[localI],        // current point label
                            savedPointLabels_[localI]   // point label when it
                                                        // was stored
                        );
                    }
                    f[fp] = addedPoints[localI];
                }
            }
        }


        // Restore
        // ~~~~~~~

        label own = mesh_.faceOwner()[masterFacei];
        label zoneID = mesh_.faceZones().whichZone(masterFacei);
        bool zoneFlip = false;
        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];
            zoneFlip = fZone.flipMap()[fZone.whichFace(masterFacei)];
        }
        label patchi = mesh_.boundaryMesh().whichPatch(masterFacei);

        if (mesh_.boundaryMesh()[patchi].coupled())
        {
            FatalErrorInFunction
                << "Master face " << masterFacei << " is on coupled patch "
                << mesh_.boundaryMesh()[patchi].name()
                << abort(FatalError);
        }

        //Pout<< "Restoring new master face " << masterFacei
        //    << " to vertices " << faces[0] << endl;

        // Modify the master face.
        meshMod.setAction
        (
            polyModifyFace
            (
                faces[0],                       // original face
                masterFacei,                    // label of face
                own,                            // owner
                -1,                             // neighbour
                false,                          // face flip
                patchi,                         // patch for face
                false,                          // remove from zone
                zoneID,                         // zone for face
                zoneFlip                        // face flip in zone
            )
        );
        restoredFaces.insert(masterFacei, masterFacei);

        // Add the previously removed faces
        for (label i = 1; i < faces.size(); i++)
        {
            //Pout<< "Restoring removed face with vertices " << faces[i]
            //    << endl;

            label facei = meshMod.setAction
            (
                polyAddFace
                (
                    faces[i],        // vertices
                    own,                    // owner,
                    -1,                     // neighbour,
                    -1,                     // masterPointID,
                    -1,                     // masterEdgeID,
                    masterFacei,             // masterFaceID,
                    false,                  // flipFaceFlux,
                    patchi,                 // patchID,
                    zoneID,                 // zoneID,
                    zoneFlip                // zoneFlip
                )
            );
            restoredFaces.insert(facei, masterFacei);
        }

        // Clear out restored set
        faceSetsVertices_[setI].clear();
        masterFace_[setI] = -1;
    }
}


// ************************************************************************* //
