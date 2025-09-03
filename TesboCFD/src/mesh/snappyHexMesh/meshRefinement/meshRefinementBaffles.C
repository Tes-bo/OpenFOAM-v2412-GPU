/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "meshRefinement.H"
#include "fvMesh.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "faceSet.H"
#include "indirectPrimitivePatch.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyModifyCell.H"
#include "polyAddFace.H"
#include "polyRemoveFace.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "unitConversion.H"
#include "OBJstream.H"
#include "patchFaceOrientation.H"
#include "PatchEdgeFaceWave.H"
#include "edgeTopoDistanceData.H"
#include "polyMeshAdder.H"
#include "IOmanip.H"
#include "refinementParameters.H"
#include "shellSurfaces.H"
#include "zeroGradientFvPatchFields.H"
#include "volFields.H"
#include "holeToFace.H"

#include "FaceCellWave.H"
#include "wallPoints.H"
#include "searchableSurfaces.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshRefinement::createBaffle
(
    const label faceI,
    const label ownPatch,
    const label neiPatch,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[faceI];
    label zoneID = mesh_.faceZones().whichZone(faceI);
    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            f,                          // modified face
            faceI,                      // label of face
            mesh_.faceOwner()[faceI],   // owner
            -1,                         // neighbour
            false,                      // face flip
            ownPatch,                   // patch for face
            false,                      // remove from zone
            zoneID,                     // zone for face
            zoneFlip                    // face flip in zone
        )
    );


    label dupFaceI = -1;

    if (mesh_.isInternalFace(faceI))
    {
        if (neiPatch == -1)
        {
            FatalErrorInFunction
                << "No neighbour patch for internal face " << faceI
                << " fc:" << mesh_.faceCentres()[faceI]
                << " ownPatch:" << ownPatch << abort(FatalError);
        }

        bool reverseFlip = false;
        if (zoneID >= 0)
        {
            reverseFlip = !zoneFlip;
        }

        dupFaceI = meshMod.setAction
        (
            polyAddFace
            (
                f.reverseFace(),            // modified face
                mesh_.faceNeighbour()[faceI],// owner
                -1,                         // neighbour
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceI,                      // masterFaceID,
                true,                       // face flip
                neiPatch,                   // patch for face
                zoneID,                     // zone for face
                reverseFlip                 // face flip in zone
            )
        );
    }
    return dupFaceI;
}


//// Check if we are a boundary face and normal of surface does
//// not align with test vector. In this case there'd probably be
//// a freestanding 'baffle' so we might as well not create it.
//// Note that since it is not a proper baffle we cannot detect it
//// afterwards so this code cannot be merged with the
//// filterDuplicateFaces code.
//bool Foam::meshRefinement::validBaffleTopology
//(
//    const label faceI,
//    const vector& n1,
//    const vector& n2,
//    const vector& testDir
//) const
//{
//
//    label patchI = mesh_.boundaryMesh().whichPatch(faceI);
//    if (patchI == -1 || mesh_.boundaryMesh()[patchI].coupled())
//    {
//        return true;
//    }
//    else if (mag(n1&n2) > cos(degToRad(30.0)))
//    {
//        // Both normals aligned. Check that test vector perpendicularish to
//        // surface normal
//        scalar magTestDir = mag(testDir);
//        if (magTestDir > VSMALL)
//        {
//            if (mag(n1&(testDir/magTestDir)) < cos(degToRad(45.0)))
//            {
//                //Pout<< "** disabling baffling face "
//                //    << mesh_.faceCentres()[faceI] << endl;
//                return false;
//            }
//        }
//    }
//    return true;
//}


void Foam::meshRefinement::getIntersections
(
    const labelList& surfacesToTest,
    const pointField& neiCc,
    const labelList& testFaces,

    labelList& globalRegion1,
    labelList& globalRegion2
) const
{
    autoPtr<OBJstream> str;
    if (debug&OBJINTERSECTIONS)
    {
        mkDir(mesh_.time().path()/timeName());
        str.reset
        (
            new OBJstream
            (
                mesh_.time().path()/timeName()/"intersections.obj"
            )
        );

        Pout<< "getIntersections : Writing surface intersections to file "
            << str().name() << nl << endl;
    }


    globalRegion1.setSize(mesh_.nFaces());
    globalRegion1 = -1;
    globalRegion2.setSize(mesh_.nFaces());
    globalRegion2 = -1;


    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    {
        labelList minLevel;
        calcCellCellRays
        (
            neiCc,
            labelList(neiCc.size(), -1),
            testFaces,
            start,
            end,
            minLevel
        );
    }


    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    surfaces_.findNearestIntersection
    (
        surfacesToTest,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
    );


    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        if (hit1[i].hit() && hit2[i].hit())
        {
            if (str)
            {
                str().writeLine(start[i], hit1[i].point());
                str().writeLine(hit1[i].point(), hit2[i].point());
                str().writeLine(hit2[i].point(), end[i]);
            }

            // Pick up the patches
            globalRegion1[faceI] =
                surfaces_.globalRegion(surface1[i], region1[i]);
            globalRegion2[faceI] =
                surfaces_.globalRegion(surface2[i], region2[i]);

            if (globalRegion1[faceI] == -1 || globalRegion2[faceI] == -1)
            {
                FatalErrorInFunction
                    << "problem." << abort(FatalError);
            }
        }
    }
}


void Foam::meshRefinement::getBafflePatches
(
    const label nErodeCellZones,
    const labelList& globalToMasterPatch,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh,
    const bool exitIfLeakPath,
    const refPtr<coordSetWriter>& leakPathFormatter,
    refPtr<surfaceWriter>& surfFormatter,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& ownPatch,
    labelList& neiPatch
) const
{
    // This determines the patches for the intersected faces to
    // - remove the outside of the mesh
    // - introduce baffles for (non-faceZone) intersections
    // Any baffles for faceZones (faceType 'baffle'/'boundary') get introduced
    // later


    // 1. Determine intersections with unnamed surfaces and cell zones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Notice that this also does hole-closure so the unnamed* is not just
    // the surface intersections.

    labelList cellToZone;
    labelList unnamedRegion1;
    labelList unnamedRegion2;
    labelList namedSurfaceRegion;
    {
        bitSet posOrientation;
        zonify
        (
            true,               // allowFreeStandingZoneFaces
            nErodeCellZones,
            -2,                 // zone to put unreached cells into
            locationsInMesh,
            zonesInMesh,
            locationsOutsideMesh,
            exitIfLeakPath,
            leakPathFormatter,
            surfFormatter,

            cellToZone,
            unnamedRegion1,
            unnamedRegion2,
            namedSurfaceRegion,
            posOrientation
        );
    }


    // 2. Baffle all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // The logic is that all intersections with unnamed surfaces become
    // baffles except where we're inbetween a cellZone and background
    // or inbetween two different cellZones.

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    ownPatch.setSize(mesh_.nFaces());
    ownPatch = -1;
    neiPatch.setSize(mesh_.nFaces());
    neiPatch = -1;

    forAll(ownPatch, faceI)
    {
        if (unnamedRegion1[faceI] != -1 || unnamedRegion2[faceI] != -1)
        {
            label ownMasterPatch = -1;
            if (unnamedRegion1[faceI] != -1)
            {
                ownMasterPatch = globalToMasterPatch[unnamedRegion1[faceI]];
            }
            label neiMasterPatch = -1;
            if (unnamedRegion2[faceI] != -1)
            {
                neiMasterPatch = globalToMasterPatch[unnamedRegion2[faceI]];
            }


            // Now we always want to produce a baffle except if
            // one side is a valid cellZone

            label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
            label neiZone = -1;

            if (mesh_.isInternalFace(faceI))
            {
                neiZone = cellToZone[mesh_.faceNeighbour()[faceI]];
            }
            else
            {
                neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];
            }


            if
            (
                (ownZone != neiZone)
             && (
                    (ownZone >= 0 && neiZone != -2)
                 || (neiZone >= 0 && ownZone != -2)
                )
             && (
                    namedSurfaceRegion.size() == 0
                 || namedSurfaceRegion[faceI] == -1
                )
            )
            {
                // One side is a valid cellZone and the neighbour is a different
                // one (or -1 but not -2). Do not baffle unless the user has
                // put both an unnamed and a named surface there. In that
                // case assume that the user wants both: baffle and faceZone.
            }
            else
            {
                ownPatch[faceI] = ownMasterPatch;
                neiPatch[faceI] = neiMasterPatch;
            }
        }
    }

    // No need to parallel sync since intersection data (surfaceIndex_ etc.)
    // already guaranteed to be synced...
    // However:
    // - owncc and neicc are reversed on different procs so might pick
    //   up different regions reversed? No problem. Neighbour on one processor
    //   might not be owner on the other processor but the neighbour is
    //   not used when creating baffles from proc faces.
    // - tolerances issues occasionally crop up.
    syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    syncTools::syncFaceList(mesh_, neiPatch, maxEqOp<label>());
}


Foam::Map<Foam::labelPair> Foam::meshRefinement::getZoneBafflePatches
(
    const bool allowBoundary,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
) const
{
    Map<labelPair> bafflePatch(mesh_.nFaces()/1000);

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    const faceZoneMesh& fZones = mesh_.faceZones();

    forAll(surfZones, surfI)
    {
        const wordList& faceZoneNames = surfZones[surfI].faceZoneNames();

        forAll(faceZoneNames, fzi)
        {
            // Get zone
            const word& faceZoneName = faceZoneNames[fzi];
            label zoneI = fZones.findZoneID(faceZoneName);
            const faceZone& fZone = fZones[zoneI];

            // Get patch allocated for zone
            label globalRegionI = surfaces_.globalRegion(surfI, fzi);
            labelPair zPatches
            (
                globalToMasterPatch[globalRegionI],
                globalToSlavePatch[globalRegionI]
            );

            Info<< "For zone " << fZone.name() << " found patches "
                << mesh_.boundaryMesh()[zPatches[0]].name() << " and "
                << mesh_.boundaryMesh()[zPatches[1]].name()
                << endl;

            forAll(fZone, i)
            {
                label faceI = fZone[i];

                if (allowBoundary || mesh_.isInternalFace(faceI))
                {
                    labelPair patches = zPatches;
                    if (fZone.flipMap()[i])
                    {
                       patches = reverse(patches);
                    }

                    if (!bafflePatch.insert(faceI, patches))
                    {
                        FatalErrorInFunction
                            << "Face " << faceI
                            << " fc:" << mesh_.faceCentres()[faceI]
                            << " in zone " << fZone.name()
                            << " is in multiple zones!"
                            << abort(FatalError);
                    }
                }
            }
        }
    }
    return bafflePatch;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createBaffles
(
    const labelList& ownPatch,
    const labelList& neiPatch
)
{
    if
    (
        ownPatch.size() != mesh_.nFaces()
     || neiPatch.size() != mesh_.nFaces()
    )
    {
        FatalErrorInFunction
            << "Illegal size :"
            << " ownPatch:" << ownPatch.size()
            << " neiPatch:" << neiPatch.size()
            << ". Should be number of faces:" << mesh_.nFaces()
            << abort(FatalError);
    }

    if (debug)
    {
        labelList syncedOwnPatch(ownPatch);
        syncTools::syncFaceList(mesh_, syncedOwnPatch, maxEqOp<label>());
        labelList syncedNeiPatch(neiPatch);
        syncTools::syncFaceList(mesh_, syncedNeiPatch, maxEqOp<label>());

        forAll(syncedOwnPatch, faceI)
        {
            if
            (
                (ownPatch[faceI] == -1 && syncedOwnPatch[faceI] != -1)
             || (neiPatch[faceI] == -1 && syncedNeiPatch[faceI] != -1)
            )
            {
                FatalErrorInFunction
                    << "Non synchronised at face:" << faceI
                    << " on patch:" << mesh_.boundaryMesh().whichPatch(faceI)
                    << " fc:" << mesh_.faceCentres()[faceI] << endl
                    << "ownPatch:" << ownPatch[faceI]
                    << " syncedOwnPatch:" << syncedOwnPatch[faceI]
                    << " neiPatch:" << neiPatch[faceI]
                    << " syncedNeiPatch:" << syncedNeiPatch[faceI]
                    << abort(FatalError);
            }
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nBaffles = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (ownPatch[faceI] != -1)
        {
            // Create baffle or repatch face. Return label of inserted baffle
            // face.
            createBaffle
            (
                faceI,
                ownPatch[faceI],   // owner side patch
                neiPatch[faceI],   // neighbour side patch
                meshMod
            );
            nBaffles++;
        }
    }
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        label coupledPatchI = -1;
        if
        (
            pp.coupled()
        && !refCast<const coupledPolyPatch>(pp).separated()
        )
        {
            coupledPatchI = patchI;
        }

        forAll(pp, i)
        {
            label faceI = pp.start()+i;

            if (ownPatch[faceI] != -1)
            {
                createBaffle
                (
                    faceI,
                    ownPatch[faceI],    // owner side patch
                    neiPatch[faceI],    // neighbour side patch
                    meshMod
                );

                if (coupledPatchI != -1)
                {
                    faceToCoupledPatch_.insert(faceI, coupledPatchI);
                }

                nBaffles++;
            }
        }
    }


    autoPtr<mapPolyMesh> mapPtr;
    if (returnReduceOr(nBaffles))
    {
        // Remove any unnecessary fields
        mesh_.clearOut();
        mesh_.moving(false);

        // Change the mesh (no inflation, parallel sync)
        mapPtr = meshMod.changeMesh(mesh_, false, true);
        mapPolyMesh& map = *mapPtr;

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map.hasMotionPoints())
        {
            mesh_.movePoints(map.preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }


        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        //- Redo the intersections on the newly create baffle faces. Note that
        //  this changes also the cell centre positions.
        faceSet baffledFacesSet(mesh_, "baffledFacesSet", 2*nBaffles);

        const labelList& reverseFaceMap = map.reverseFaceMap();
        const labelList& faceMap = map.faceMap();

        // Pick up owner side of baffle
        forAll(ownPatch, oldFaceI)
        {
            label faceI = reverseFaceMap[oldFaceI];

            if (ownPatch[oldFaceI] != -1 && faceI >= 0)
            {
                const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

                baffledFacesSet.insert(ownFaces);
            }
        }
        // Pick up neighbour side of baffle (added faces)
        forAll(faceMap, faceI)
        {
            label oldFaceI = faceMap[faceI];

            if (oldFaceI >= 0 && reverseFaceMap[oldFaceI] != faceI)
            {
                const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

                baffledFacesSet.insert(ownFaces);
            }
        }
        baffledFacesSet.sync(mesh_);

        updateMesh(map, baffledFacesSet.toc());
    }

    return mapPtr;
}


Foam::labelList Foam::meshRefinement::getZones
(
    const List<surfaceZonesInfo::faceZoneType>& fzTypes
) const
{
    const faceZoneMesh& faceZones = mesh_.faceZones();

    DynamicList<label> zoneIDs(faceZones.size());

    forAll(faceZones, zoneI)
    {
        const faceZone& fZone = faceZones[zoneI];

        label mpI, spI;
        surfaceZonesInfo::faceZoneType fzType;
        bool hasInfo = getFaceZoneInfo(fZone.name(), mpI, spI, fzType);

        if (hasInfo && fzTypes.found(fzType))
        {
            zoneIDs.append(zoneI);
        }
    }
    return zoneIDs;
}


// Subset those baffles where both faces are on the same zone
Foam::List<Foam::labelPair> Foam::meshRefinement::subsetBaffles
(
    const polyMesh& mesh,
    const labelList& zoneIDs,
    const List<labelPair>& baffles
)
{
    const faceZoneMesh& faceZones = mesh.faceZones();

    // Mark zone per face
    labelList faceToZone(mesh.nFaces(), -1);

    for (const label zoneID : zoneIDs)
    {
        labelUIndList(faceToZone, faceZones[zoneID]) = zoneID;
    }


    // Subset baffles
    DynamicList<labelPair> newBaffles(baffles.size());
    forAll(baffles, i)
    {
        const labelPair& p = baffles[i];
        if (faceToZone[p[0]] != -1 && (faceToZone[p[0]] == faceToZone[p[1]]))
        {
            newBaffles.append(p);
        }
    }

    return newBaffles;
}


void Foam::meshRefinement::mapBaffles
(
    const polyMesh& mesh,
    const labelList& faceMap,
    List<labelPair>& baffles
)
{
    // Create old-to-new map just for boundary faces. (since multiple faces
    // get created from the same baffle face)
    labelList reverseFaceMap(mesh.nFaces(), -1);
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        label oldFacei = faceMap[facei];
        if (oldFacei != -1)
        {
            reverseFaceMap[oldFacei] = facei;
        }
    }


    DynamicList<labelPair> newBaffles(baffles.size());
    forAll(baffles, i)
    {
        const labelPair& p = baffles[i];
        labelPair newBaffle
        (
            reverseFaceMap[p[0]],
            reverseFaceMap[p[1]]
        );
        if (newBaffle[0] != -1 && newBaffle[1] != -1)
        {
            newBaffles.append(newBaffle);
        }
    }
    baffles = std::move(newBaffles);
}


void Foam::meshRefinement::getZoneFaces
(
    const labelList& zoneIDs,
    labelList& faceZoneID,
    labelList& ownPatch,
    labelList& neiPatch,
    labelList& nBaffles
) const
{
    const faceZoneMesh& faceZones = mesh_.faceZones();

    // Per (internal) face the patch related to the faceZone
    ownPatch.setSize(mesh_.nFaces());
    ownPatch= -1;
    neiPatch.setSize(mesh_.nFaces());
    neiPatch = -1;
    faceZoneID.setSize(mesh_.nFaces());
    faceZoneID = -1;
    nBaffles.setSize(zoneIDs.size());
    nBaffles = Zero;


    //- Get per face whether it is internal or coupled
    const bitSet isInternalOrCoupled
    (
        syncTools::getInternalOrCoupledFaces(mesh_)
    );

    forAll(zoneIDs, j)
    {
        label zoneI = zoneIDs[j];
        const faceZone& fz = faceZones[zoneI];
        const word& masterName = faceZoneToMasterPatch_[fz.name()];
        label masterPatchI = mesh_.boundaryMesh().findPatchID(masterName);
        const word& slaveName = faceZoneToSlavePatch_[fz.name()];
        label slavePatchI = mesh_.boundaryMesh().findPatchID(slaveName);

        if (masterPatchI == -1 || slavePatchI == -1)
        {
            FatalErrorInFunction
                << "Problem: masterPatchI:" << masterPatchI
                << " slavePatchI:" << slavePatchI << exit(FatalError);
        }

        forAll(fz, i)
        {
            label faceI = fz[i];
            if (isInternalOrCoupled[faceI])
            {
                if (fz.flipMap()[i])
                {
                    ownPatch[faceI] = slavePatchI;
                    neiPatch[faceI] = masterPatchI;
                }
                else
                {
                    ownPatch[faceI] = masterPatchI;
                    neiPatch[faceI] = slavePatchI;
                }
                faceZoneID[faceI] = zoneI;

                nBaffles[j]++;
            }
        }
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createZoneBaffles
(
    const labelList& zoneIDs,
    List<labelPair>& baffles,
    labelList& originatingFaceZone
)
{
    autoPtr<mapPolyMesh> map;

    if (zoneIDs.size() > 0)
    {
        const faceZoneMesh& faceZones = mesh_.faceZones();

        // Split internal faces on interface surfaces
        Info<< "Converting zoned faces into baffles ..." << endl;

        // Get faceZone and patch(es) per face (or -1 if face not on faceZone)
        labelList faceZoneID;
        labelList ownPatch;
        labelList neiPatch;
        labelList nBaffles;
        getZoneFaces(zoneIDs, faceZoneID, ownPatch, neiPatch, nBaffles);

        label nLocalBaffles = sum(nBaffles);

        label nTotalBaffles = returnReduce(nLocalBaffles, sumOp<label>());

        if (nTotalBaffles > 0)
        {
            Pstream::listCombineReduce(nBaffles, plusEqOp<label>());

            Info<< nl
                << setf(ios_base::left)
                << setw(30) << "FaceZone"
                << setw(10) << "FaceType"
                << setw(10) << "nBaffles"
                << nl
                << setw(30) << "--------"
                << setw(10) << "--------"
                << setw(10) << "--------"
                << endl;

            forAll(zoneIDs, j)
            {
                label zoneI = zoneIDs[j];
                const faceZone& fz = faceZones[zoneI];

                label mpI, spI;
                surfaceZonesInfo::faceZoneType fzType;
                bool hasInfo = getFaceZoneInfo(fz.name(), mpI, spI, fzType);

                if (hasInfo)
                {
                    Info<< setf(ios_base::left)
                        << setw(30) << fz.name()
                        << setw(10)
                        << surfaceZonesInfo::faceZoneTypeNames[fzType]
                        << setw(10) << nBaffles[j]
                        << nl;
                }
            }
            Info<< endl;

            // Create baffles.
            map = createBaffles(ownPatch, neiPatch);

            // Get pairs of faces created.
            // Just loop over faceMap and store baffle if we encounter a slave
            // face.

            baffles.setSize(nLocalBaffles);
            originatingFaceZone.setSize(nLocalBaffles);
            label baffleI = 0;

            const labelList& faceMap = map().faceMap();
            const labelList& reverseFaceMap = map().reverseFaceMap();

            for
            (
                label faceI = mesh_.nInternalFaces();
                faceI < mesh_.nFaces();
                faceI++
            )
            {
                label oldFaceI = faceMap[faceI];
                label masterFaceI = reverseFaceMap[oldFaceI];
                if (masterFaceI != faceI && ownPatch[oldFaceI] != -1)
                {
                    baffles[baffleI] = labelPair(masterFaceI, faceI);
                    originatingFaceZone[baffleI] = faceZoneID[oldFaceI];
                    baffleI++;
                }
            }

            if (baffleI != baffles.size())
            {
                FatalErrorInFunction
                    << "Had " << baffles.size() << " baffles to create "
                    << " but encountered " << baffleI
                    << " slave faces originating from patcheable faces."
                    << abort(FatalError);
            }

            if (debug&MESH)
            {
                const_cast<Time&>(mesh_.time())++;
                Pout<< "Writing zone-baffled mesh to time " << timeName()
                    << endl;
                write
                (
                    debugType(debug),
                    writeType(writeLevel() | WRITEMESH),
                    mesh_.time().path()/"baffles"
                );
            }
        }
        Info<< "Created " << nTotalBaffles << " baffles in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
    else
    {
        baffles.clear();
        originatingFaceZone.clear();
    }

    return map;
}


Foam::List<Foam::labelPair> Foam::meshRefinement::freeStandingBaffles
(
    const List<labelPair>& couples,
    const scalar planarAngle,
    const bool samePatch
) const
{
    // Done by counting the number of baffles faces per mesh edge. If edge
    // has 2 boundary faces and both are baffle faces it is the edge of a baffle
    // region.

    // All duplicate faces on edge of the patch are to be merged.
    // So we count for all edges of duplicate faces how many duplicate
    // faces use them.
    labelList nBafflesPerEdge(mesh_.nEdges(), Zero);


    // This algorithm is quite tricky. We don't want to use edgeFaces and
    // also want it to run in parallel so it is now an algorithm over
    // all (boundary) faces instead.
    // We want to pick up any edges that are only used by the baffle
    // or internal faces but not by any other boundary faces. So
    // - increment count on an edge by 1 if it is used by any (uncoupled)
    //   boundary face.
    // - increment count on an edge by 1000000 if it is used by a baffle face
    // - sum in parallel
    //
    // So now any edge that is used by baffle faces only will have the
    // value 2*1000000+2*1.


    const label baffleValue = 1000000;


    // Count number of boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Count number of boundary faces. Discard coupled boundary faces.
        if (!pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                const labelList& fEdges = mesh_.faceEdges(faceI);

                forAll(fEdges, fEdgeI)
                {
                    nBafflesPerEdge[fEdges[fEdgeI]]++;
                }
                faceI++;
            }
        }
    }


    DynamicList<label> fe0;
    DynamicList<label> fe1;


    // Count number of duplicate boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(couples, i)
    {
        {
            label f0 = couples[i].first();
            const labelList& fEdges0 = mesh_.faceEdges(f0, fe0);
            forAll(fEdges0, fEdgeI)
            {
                nBafflesPerEdge[fEdges0[fEdgeI]] += baffleValue;
            }
        }

        {
            label f1 = couples[i].second();
            const labelList& fEdges1 = mesh_.faceEdges(f1, fe1);
            forAll(fEdges1, fEdgeI)
            {
                nBafflesPerEdge[fEdges1[fEdgeI]] += baffleValue;
            }
        }
    }

    // Add nBaffles on shared edges
    syncTools::syncEdgeList
    (
        mesh_,
        nBafflesPerEdge,
        plusEqOp<label>(),  // in-place add
        label(0)            // initial value
    );


    // Baffles which are not next to other boundaries and baffles will have
    // nBafflesPerEdge value 2*baffleValue+2*1 (from 2 boundary faces which
    // are both baffle faces)

    List<labelPair> filteredCouples(couples.size());
    label filterI = 0;

    forAll(couples, i)
    {
        const labelPair& couple = couples[i];

        if
        (
           !samePatch
         || (
                patches.whichPatch(couple.first())
             == patches.whichPatch(couple.second())
            )
        )
        {
            const labelList& fEdges = mesh_.faceEdges(couple.first());

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (nBafflesPerEdge[edgeI] == 2*baffleValue+2*1)
                {
                    filteredCouples[filterI++] = couple;
                    break;
                }
            }
        }
    }
    filteredCouples.setSize(filterI);


    label nFiltered = returnReduce(filteredCouples.size(), sumOp<label>());

    Info<< "freeStandingBaffles : detected "
        << nFiltered
        << " free-standing baffles out of "
        << returnReduce(couples.size(), sumOp<label>())
        << nl << endl;


    if (nFiltered > 0)
    {
        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(filteredCouples.size());
        pointField end(filteredCouples.size());

        const pointField& cellCentres = mesh_.cellCentres();

        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];
            start[i] = cellCentres[mesh_.faceOwner()[couple.first()]];
            end[i] = cellCentres[mesh_.faceOwner()[couple.second()]];
        }

        // Extend segments a bit
        {
            const vectorField smallVec(ROOTSMALL*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do test for intersections
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;

        surfaces_.findNearestIntersection
        (
            identity(surfaces_.surfaces().size()),
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );

        //mkDir(mesh_.time().path()/timeName());
        //OBJstream str
        //(
        //    mesh_.time().path()/timeName()/"flatBaffles.obj"
        //);

        const scalar planarAngleCos = Foam::cos(degToRad(planarAngle));

        label filterI = 0;
        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];

            // Note: for a baffle-surface we do not want to merge the baffle.
            // We could either check for hitting the same triangle (but you
            // might hit same point on neighbouring triangles due to tolerance
            // issues) or better just to compare the hit point.
            // This might still go wrong for a ray in the plane of the triangle
            // which might hit two different points on the same triangle due
            // to tolerances...

            if
            (
                hit1[i].hit()
             && hit2[i].hit()
             && hit1[i].point().dist(hit2[i].point()) > mergeDistance_
            )
            {
                // Two different hits. Check angle.
                //str.write
                //(
                //    linePointRef(hit1[i].point(), hit2[i].point()),
                //    normal1[i],
                //    normal2[i]
                //);

                if ((normal1[i]&normal2[i]) > planarAngleCos)
                {
                    // Both normals aligned
                    vector n = end[i]-start[i];
                    scalar magN = mag(n);
                    if (magN > VSMALL)
                    {
                        filteredCouples[filterI++] = couple;
                    }
                }
            }
            else if (hit1[i].hit() || hit2[i].hit())
            {
                // Single hit. Do not include in freestanding baffles.
            }
        }

        filteredCouples.setSize(filterI);

        Info<< "freeStandingBaffles : detected "
            << returnReduce(filterI, sumOp<label>())
            << " planar (within " << planarAngle
            << " degrees) free-standing baffles out of "
            << nFiltered
            << nl << endl;
    }

    return filteredCouples;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeBaffles
(
    const List<labelPair>& couples,
    const Map<label>& faceToPatch
)
{
    autoPtr<mapPolyMesh> mapPtr;

    if (returnReduceOr(couples.size() || faceToPatch.size()))
    {
        // Mesh change engine
        polyTopoChange meshMod(mesh_);

        const faceList& faces = mesh_.faces();
        const labelList& faceOwner = mesh_.faceOwner();
        const faceZoneMesh& faceZones = mesh_.faceZones();

        forAll(couples, i)
        {
            label face0 = couples[i].first();
            label face1 = couples[i].second();

            // face1 < 0 signals a coupled face that has been converted to
            // baffle

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (face1 < 0 || own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = faceZones.whichZone(face0);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                }

                label nei = (face1 < 0 ? -1 : own1);

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        nei,                    // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
            else
            {
                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(face1);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                }

                meshMod.setAction(polyRemoveFace(face0));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face1],           // modified face
                        face1,                  // label of face being modified
                        own1,                   // owner
                        own0,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }

        forAllConstIters(faceToPatch, iter)
        {
            const label faceI = iter.key();
            const label patchI = iter.val();

            if (!mesh_.isInternalFace(faceI))
            {
                FatalErrorInFunction
                    << "problem: face:" << faceI
                    << " at:" << mesh_.faceCentres()[faceI]
                    << "(wanted patch:" << patchI
                    << ") is an internal face" << exit(FatalError);
            }

            label zoneID = faceZones.whichZone(faceI);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    faces[faceI],           // modified face
                    faceI,                  // label of face being modified
                    faceOwner[faceI],       // owner
                    -1,                     // neighbour
                    false,                  // face flip
                    patchI,                 // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }


        // Remove any unnecessary fields
        mesh_.clearOut();
        mesh_.moving(false);

        // Change the mesh (no inflation)
        mapPtr = meshMod.changeMesh(mesh_, false, true);
        mapPolyMesh& map = *mapPtr;

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map.hasMotionPoints())
        {
            mesh_.movePoints(map.preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Recalculate intersections on merged faces since
        // this seems to give problems? Note: should not be necessary since
        // baffles preserve intersections from when they were created.
        labelList newExposedFaces(2*couples.size());
        label newI = 0;

        forAll(couples, i)
        {
            const label newFace0 = map.reverseFaceMap()[couples[i].first()];
            if (newFace0 != -1)
            {
                newExposedFaces[newI++] = newFace0;
            }

            const label newFace1 = map.reverseFaceMap()[couples[i].second()];
            if (newFace1 != -1)
            {
                newExposedFaces[newI++] = newFace1;
            }
        }
        newExposedFaces.setSize(newI);
        updateMesh(map, newExposedFaces);
    }

    return mapPtr;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeZoneBaffles
(
    const bool doInternalZones,
    const bool doBaffleZones
)
{
    labelList zoneIDs;
    {
        DynamicList<surfaceZonesInfo::faceZoneType> fzTypes;
        if (doInternalZones)
        {
            fzTypes.append(surfaceZonesInfo::INTERNAL);
        }
        if (doBaffleZones)
        {
            fzTypes.append(surfaceZonesInfo::BAFFLE);
        }
        zoneIDs = getZones(fzTypes);
    }

    List<labelPair> zoneBaffles
    (
        subsetBaffles
        (
            mesh_,
            zoneIDs,
            localPointRegion::findDuplicateFacePairs(mesh_)
        )
    );

    autoPtr<mapPolyMesh> mapPtr;
    if (returnReduceOr(zoneBaffles.size()))
    {
        mapPtr = mergeBaffles(zoneBaffles, Map<label>(0));
    }
    return mapPtr;
}


// Finds region per cell for cells inside closed named surfaces
void Foam::meshRefinement::findCellZoneGeometric
(
    const pointField& neiCc,
    const labelList& closedNamedSurfaces,   // indices of closed surfaces
    labelList& namedSurfaceRegion,          // per face: named surface region
    const labelList& surfaceToCellZone,     // cell zone index per surface

    labelList& cellToZone
) const
{
    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Check if cell centre is inside
    labelList insideSurfaces;
    surfaces_.findInside
    (
        closedNamedSurfaces,
        cellCentres,
        insideSurfaces
    );

    forAll(insideSurfaces, cellI)
    {
        label surfI = insideSurfaces[cellI];

        if (surfI != -1)
        {
            if (cellToZone[cellI] == -2)
            {
                cellToZone[cellI] = surfaceToCellZone[surfI];
            }
            else if (cellToZone[cellI] == -1)
            {
                // ? Allow named surface to override background zone (-1)
                // This is used in the multiRegionHeater tutorial where the
                // locationInMesh is inside a named surface.
                cellToZone[cellI] = surfaceToCellZone[surfI];
            }
        }
    }


    // Some cells with cell centres close to surface might have
    // had been put into wrong surface. Recheck with perturbed cell centre.


    // 1. Collect points

    // Count points to test.
    label nCandidates = 0;
    forAll(namedSurfaceRegion, faceI)
    {
        if (namedSurfaceRegion[faceI] != -1)
        {
            if (mesh_.isInternalFace(faceI))
            {
                nCandidates += 2;
            }
            else
            {
                nCandidates += 1;
            }
        }
    }

    // Collect points.
    pointField candidatePoints(nCandidates);
    nCandidates = 0;
    forAll(namedSurfaceRegion, faceI)
    {
        if (namedSurfaceRegion[faceI] != -1)
        {
            label own = faceOwner[faceI];
            const point& ownCc = cellCentres[own];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = faceNeighbour[faceI];
                const point& neiCc = cellCentres[nei];
                // Perturbed cc
                const vector d = 1e-4*(neiCc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
                candidatePoints[nCandidates++] = neiCc+d;
            }
            else
            {
                //const point& neiFc = mesh_.faceCentres()[faceI];
                const point& neiFc = neiCc[faceI-mesh_.nInternalFaces()];

                // Perturbed cc
                const vector d = 1e-4*(neiFc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
            }
        }
    }


    // 2. Test points for inside

    surfaces_.findInside
    (
        closedNamedSurfaces,
        candidatePoints,
        insideSurfaces
    );


    // 3. Update zone information

    nCandidates = 0;
    forAll(namedSurfaceRegion, faceI)
    {
        if (namedSurfaceRegion[faceI] != -1)
        {
            label own = faceOwner[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }

                label neiSurfI = insideSurfaces[nCandidates++];
                if (neiSurfI != -1)
                {
                    label nei = faceNeighbour[faceI];

                    cellToZone[nei] = surfaceToCellZone[neiSurfI];
                }
            }
            else
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }
            }
        }
    }


    // Adapt the namedSurfaceRegion
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // for if any cells were not completely covered.

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
        label neiZone = cellToZone[mesh_.faceNeighbour()[faceI]];

        if (namedSurfaceRegion[faceI] == -1 && (ownZone != neiZone))
        {
            // Give face the zone of min cell zone (but only if the
            // cellZone originated from a closed, named surface)

            label minZone;
            if (ownZone == -1)
            {
                minZone = neiZone;
            }
            else if (neiZone == -1)
            {
                minZone = ownZone;
            }
            else
            {
                minZone = min(ownZone, neiZone);
            }

            // Make sure the cellZone originated from a closed surface. Use
            // hardcoded region 0 inside named surface.
            label geomSurfI = surfaceToCellZone.find(minZone);

            if (geomSurfI != -1)
            {
                namedSurfaceRegion[faceI] =
                    surfaces_.globalRegion(geomSurfI, 0);
            }
        }
    }

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];

                if (namedSurfaceRegion[faceI] == -1 && (ownZone != neiZone))
                {
                    // Give face the min cell zone
                    label minZone;
                    if (ownZone == -1)
                    {
                        minZone = neiZone;
                    }
                    else if (neiZone == -1)
                    {
                        minZone = ownZone;
                    }
                    else
                    {
                        minZone = min(ownZone, neiZone);
                    }

                    // Make sure the cellZone originated from a closed surface.
                    // Use hardcoded region 0 inside named surface.
                    label geomSurfI = surfaceToCellZone.find(minZone);

                    if (geomSurfI != -1)
                    {
                        namedSurfaceRegion[faceI] =
                            surfaces_.globalRegion(geomSurfI, 0);
                    }
                }
            }
        }
    }

    // Sync
    syncTools::syncFaceList(mesh_, namedSurfaceRegion, maxEqOp<label>());
}


void Foam::meshRefinement::findCellZoneInsideWalk
(
    const pointField& locationsInMesh,
    const labelList& zonesInMesh,
    const labelList& faceToZone, // per face -1 or some index >= 0

    labelList& cellToZone
) const
{
    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());
    //selectSeparatedCoupledFaces(blockedFace);

    forAll(blockedFace, faceI)
    {
        if (faceToZone[faceI] == -1)
        {
            blockedFace[faceI] = false;
        }
        else
        {
            blockedFace[faceI] = true;
        }
    }
    // No need to sync since faceToZone already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];
        label zoneID = zonesInMesh[i];

        // Find the region containing the insidePoint
        label keepRegionI = findRegion
        (
            mesh_,
            cellRegion,
            vector::uniform(mergeDistance_),
            insidePoint
        );

        Info<< "For cellZone "
            <<  (
                    zoneID == -1
                  ? "none"
                  : mesh_.cellZones()[zoneID].name()
                )
            << " found point " << insidePoint
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }



        // Set all cells with this region to the zoneID
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        label nWarnings = 0;

        forAll(cellRegion, cellI)
        {
            if (cellRegion[cellI] == keepRegionI)
            {
                if (cellToZone[cellI] == -2)
                {
                    // First visit of cell
                    cellToZone[cellI] = zoneID;
                }
                else if (cellToZone[cellI] != zoneID)
                {
                    if (cellToZone[cellI] != -1 && nWarnings < 10)
                    {
                        WarningInFunction
                            << "Cell " << cellI
                            << " at " << mesh_.cellCentres()[cellI]
                            << " is inside cellZone "
                            <<  (
                                    zoneID == -1
                                  ? "none"
                                  : mesh_.cellZones()[zoneID].name()
                                )
                            << " from locationInMesh " << insidePoint
                            << " but already marked as being in zone "
                            << mesh_.cellZones()[cellToZone[cellI]].name()
                            << endl
                            << "This can happen if your surfaces are not"
                            << " (sufficiently) closed."
                            << endl;
                        nWarnings++;
                    }

                    // Override the zone
                    cellToZone[cellI] = zoneID;
                }
            }
        }
    }
}


void Foam::meshRefinement::findCellZoneInsideWalk
(
    const pointField& locationsInMesh,
    const wordList& zoneNamesInMesh,
    const labelList& faceToZone,    // per face -1 or some index >= 0

    labelList& cellToZone
) const
{
    const cellZoneMesh& czs = mesh_.cellZones();

    labelList zoneIDs(zoneNamesInMesh.size());
    forAll(zoneNamesInMesh, i)
    {
        zoneIDs[i] = czs.findZoneID(zoneNamesInMesh[i]);
    }
    findCellZoneInsideWalk
    (
        locationsInMesh,
        zoneIDs,
        faceToZone,
        cellToZone
    );
}


bool Foam::meshRefinement::calcRegionToZone
(
    const label backgroundZoneID,
    const label surfZoneI,
    const label ownRegion,
    const label neiRegion,

    labelList& regionToCellZone
) const
{
    bool changed = false;

    // Check whether inbetween different regions
    if (ownRegion != neiRegion)
    {
        // Jump. Change one of the sides to my type.

        // 1. Interface between my type and unset region.
        // Set region to keepRegion

        if (regionToCellZone[ownRegion] == -2)
        {
            if (surfZoneI == -1)
            {
                // Special: face is -on faceZone  -not real boundary
                //          -not on cellZone
                // so make regions same on either side
                if (regionToCellZone[neiRegion] != -2)
                {
                    regionToCellZone[ownRegion] = regionToCellZone[neiRegion];
                    changed = true;
                }
            }
            else if (regionToCellZone[neiRegion] == surfZoneI)
            {
                // Face between unset and my region. Put unset
                // region into background region
                //MEJ: see comment in findCellZoneTopo
                if (backgroundZoneID != -2)
                {
                    regionToCellZone[ownRegion] = backgroundZoneID;
                    changed = true;
                }
            }
            else if (regionToCellZone[neiRegion] != -2)
            {
                // Face between unset and other region.
                // Put unset region into my region
                regionToCellZone[ownRegion] = surfZoneI;
                changed = true;
            }
        }
        else if (regionToCellZone[neiRegion] == -2)
        {
            if (surfZoneI == -1)
            {
                // Special: face is -on faceZone  -not real boundary
                //          -not on cellZone
                // so make regions same on either side
                regionToCellZone[neiRegion] = regionToCellZone[ownRegion];
                changed = true;
            }
            else if (regionToCellZone[ownRegion] == surfZoneI)
            {
                // Face between unset and my region. Put unset
                // region into background region
                if (backgroundZoneID != -2)
                {
                    regionToCellZone[neiRegion] = backgroundZoneID;
                    changed = true;
                }
            }
            else if (regionToCellZone[ownRegion] != -2)
            {
                // Face between unset and other region.
                // Put unset region into my region
                regionToCellZone[neiRegion] = surfZoneI;
                changed = true;
            }
        }
    }
    return changed;
}


void Foam::meshRefinement::findCellZoneTopo
(
    const label backgroundZoneID,
    const pointField& locationsInMesh,
    const labelList& unnamedSurfaceRegion,
    const labelList& namedSurfaceRegion,
    const labelList& surfaceToCellZone,
    labelList& cellToZone
) const
{
    // This routine fixes small problems with left over unassigned regions
    // (after all off the unreachable bits of the mesh have been removed).
    // This routine splits the mesh into regions, based on the intersection
    // with a surface. The problem is that we know the surface which
    // (intersected) face belongs to (in namedSurfaceRegion) but we don't
    // know which side of the face it relates to. So all we are doing here
    // is get the correspondence between surface/cellZone and regionSplit
    // region. See the logic in calcRegionToZone.
    // Basically it looks at the neighbours of a face on a named surface.
    // If one side is in the cellZone corresponding to the surface and
    // the other side is unassigned (-2) it sets this to the background zone.
    // So the zone to set these unmarked cells to is provided as argument:
    // - backgroundZoneID = -2 : do not change so remove cells
    // - backgroundZoneID = -1 : put into background

    // Assumes:
    // - region containing keepPoint does not go into a cellZone
    // - all other regions can be found by crossing faces marked in
    //   namedSurfaceRegion.

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());

    forAll(unnamedSurfaceRegion, faceI)
    {
        if
        (
            unnamedSurfaceRegion[faceI] == -1
         && namedSurfaceRegion[faceI] == -1
        )
        {
            blockedFace[faceI] = false;
        }
        else
        {
            blockedFace[faceI] = true;
        }
    }
    // No need to sync since namedSurfaceRegion already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Per mesh region the zone the cell should be put in.
    // -2   : not analysed yet
    // -1   : keepPoint region. Not put into any cellzone.
    // >= 0 : index of cellZone
    labelList regionToCellZone(cellRegion.nRegions(), -2);

    // See which cells already are set in the cellToZone (from geometric
    // searching) and use these to take over their zones.
    // Note: could be improved to count number of cells per region.
    forAll(cellToZone, cellI)
    {
        label regionI = cellRegion[cellI];
        if (cellToZone[cellI] != -2 && regionToCellZone[regionI] == -2)
        {
            regionToCellZone[regionI] = cellToZone[cellI];
        }
    }

    // Synchronise regionToCellZone.
    // Note:
    // - region numbers are identical on all processors
    // - keepRegion is identical ,,
    // - cellZones are identical ,,
    Pstream::listCombineReduce(regionToCellZone, maxEqOp<label>());


    // Find the region containing the keepPoint
    forAll(locationsInMesh, i)
    {
        const point& keepPoint = locationsInMesh[i];
        label keepRegionI = findRegion
        (
            mesh_,
            cellRegion,
            vector::uniform(mergeDistance_),
            keepPoint
        );

        Info<< "Found point " << keepPoint
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << keepPoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        // Mark default region with zone -1. Note that all regions should
        // already be matched to a cellZone through the loop over cellToZone.
        // This is just to mop up any remaining regions. Not sure whether
        // this is needed?
        if (regionToCellZone[keepRegionI] == -2)
        {
            regionToCellZone[keepRegionI] = -1;
        }
    }


    // Find correspondence between cell zone and surface
    // by changing cell zone every time we cross a surface.
    while (true)
    {
        // Synchronise regionToCellZone.
        // Note:
        // - region numbers are identical on all processors
        // - keepRegion is identical ,,
        // - cellZones are identical ,,
        // This done at top of loop to account for geometric matching
        // not being synchronised.
        Pstream::listCombineReduce(regionToCellZone, maxEqOp<label>());


        bool changed = false;

        // Internal faces

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label regionI = namedSurfaceRegion[faceI];

            // Connected even if no cellZone defined for surface
            if (unnamedSurfaceRegion[faceI] == -1 && regionI != -1)
            {
                // Calculate region to zone from cellRegions on either side
                // of internal face.

                label surfI = surfaces_.whichSurface(regionI).first();

                bool changedCell = calcRegionToZone
                (
                    backgroundZoneID,
                    surfaceToCellZone[surfI],
                    cellRegion[mesh_.faceOwner()[faceI]],
                    cellRegion[mesh_.faceNeighbour()[faceI]],
                    regionToCellZone
                );

                changed = changed | changedCell;
            }
        }

        // Boundary faces

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Get coupled neighbour cellRegion
        labelList neiCellRegion;
        syncTools::swapBoundaryCellList(mesh_, cellRegion, neiCellRegion);

        // Calculate region to zone from cellRegions on either side of coupled
        // face.
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;

                    label regionI = namedSurfaceRegion[faceI];

                    // Connected even if no cellZone defined for surface
                    if (unnamedSurfaceRegion[faceI] == -1 && regionI != -1)
                    {
                        label surfI = surfaces_.whichSurface(regionI).first();

                        bool changedCell = calcRegionToZone
                        (
                            backgroundZoneID,
                            surfaceToCellZone[surfI],
                            cellRegion[mesh_.faceOwner()[faceI]],
                            neiCellRegion[faceI-mesh_.nInternalFaces()],
                            regionToCellZone
                        );

                        changed = changed | changedCell;
                    }
                }
            }
        }

        if (!returnReduceOr(changed))
        {
            break;
        }
    }


    if (debug)
    {
        Pout<< "meshRefinement::findCellZoneTopo :"
            << " nRegions:" << regionToCellZone.size()
            << " of which visited (-1 = background, >= 0 : cellZone) :"
            << endl;

        forAll(regionToCellZone, regionI)
        {
            if (regionToCellZone[regionI] != -2)
            {
                Pout<< "Region " << regionI
                    << " becomes cellZone:" << regionToCellZone[regionI]
                    << endl;
            }
        }
    }

    // Rework into cellToZone
    forAll(cellToZone, cellI)
    {
        label regionI = cellRegion[cellI];
        if (cellToZone[cellI] == -2 && regionToCellZone[regionI] != -2)
        {
            cellToZone[cellI] = regionToCellZone[regionI];
        }
    }
}


void Foam::meshRefinement::erodeCellZone
(
    const label nErodeCellZones,
    const label backgroundZoneID,
    const labelList& unnamedSurfaceRegion,
    const labelList& namedSurfaceRegion,
    labelList& cellToZone
) const
{
    // This routine fixes small problems with left over unassigned regions
    // (after all off the unreachable bits of the mesh have been removed).
    // The problem is that the cell zone information might be inconsistent
    // with the face zone information. So what we do here is to erode
    // any cell zones until we hit a named face.
    // - backgroundZoneID = -2 : do not change so remove cells
    // - backgroundZoneID = -1 : put into background
    // Note that is the opposite of findCellZoneTopo which moves unassigned
    // regions into a neighbouring region(=cellZone) unless there is an
    // intersected faces inbetween the two.

    for (label iter = 0; iter < nErodeCellZones; iter++)
    {
        label nChanged = 0;

        labelList erodedCellToZone(cellToZone);

        // Do internal faces
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if
            (
                unnamedSurfaceRegion[facei] == -1
             && namedSurfaceRegion[facei] == -1
            )
            {
                label own = mesh_.faceOwner()[facei];
                label nei = mesh_.faceNeighbour()[facei];
                if (cellToZone[own] == -2 && cellToZone[nei] >= -1)
                {
                    erodedCellToZone[nei] = backgroundZoneID;
                    nChanged++;
                }
                else if (cellToZone[nei] == -2 && cellToZone[own] >= -1)
                {
                    erodedCellToZone[own] = backgroundZoneID;
                    nChanged++;
                }
            }
        }

        // Do boundary faces

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Get coupled neighbour cellRegion
        labelList neiCellZone;
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

        // Calculate region to zone from cellRegions on either side of coupled
        // face.
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    if
                    (
                        unnamedSurfaceRegion[facei] == -1
                     && namedSurfaceRegion[facei] == -1
                    )
                    {
                        label own = mesh_.faceOwner()[facei];
                        label bFacei = facei-mesh_.nInternalFaces();
                        if (neiCellZone[bFacei] == -2 && cellToZone[own] >= -1)
                        {
                            erodedCellToZone[own] = backgroundZoneID;
                            nChanged++;
                        }
                    }
                }
            }
        }

        cellToZone.transfer(erodedCellToZone);

        reduce(nChanged, sumOp<label>());
        if (debug)
        {
            Pout<< "erodeCellZone : eroded " << nChanged
                << " cells (moved from cellZone to background zone "
                << backgroundZoneID << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }
}


void Foam::meshRefinement::growCellZone
(
    const label nGrowCellZones,
    const label backgroundZoneID,
    labelList& unnamedSurfaceRegion1,
    labelList& unnamedSurfaceRegion2,
    labelList& namedSurfaceRegion,  // potentially zero size if no faceZones
    labelList& cellToZone
) const
{
    if (nGrowCellZones != 1)
    {
        WarningInFunction
            << "Growing currently only supported for single layer."
            << " Exiting zone growing" << endl;
        return;
    }


    // See meshRefinement::markProximityRefinementWave:
    // - walk out nGrowCellZones iterations from outside of cellZone
    //   and wall into unassigned cells
    // - detect any cells inbetween
    //      - multiple zones
    //      - zone and wall
    //   and
    //      - assign cells to the cellZone
    //      - unblock faces (along the path) inbetween


    // Special tag for walls
    const FixedList<label, 3> wallTag
    ({
        labelMax,   // Special value for wall face. Loses out to cellZone tag
        labelMax,
        labelMax
    });

    // Work arrays
    pointField origins(1);
    scalarField distSqrs(1);
    List<FixedList<label, 3>> surfaces(1);


    // Set up blockage. Marked with 0 distance (so always wins)

    //List<wallPoints> allFaceInfo(mesh_.nFaces());
    //for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    //{
    //    if
    //    (
    //        unnamedSurfaceRegion1[facei] != -1
    //     || unnamedSurfaceRegion2[facei] != -1
    //    )
    //    {
    //        origins[0] = mesh_.faceCentres()[facei];
    //        distSqrs[0] = 0.0;  // Initial distance
    //        surfaces[0] = wallTag;
    //        allFaceInfo[facei] = wallPoints(origins, distSqrs, surfaces);
    //    }
    //}
    List<wallPoints> allCellInfo(mesh_.nCells());
    forAll(cellToZone, celli)
    {
        if (cellToZone[celli] >= 0)
        {
            const FixedList<label, 3> zoneTag
            ({
                cellToZone[celli],      // zone
                0,                      // 'region'
                labelMax                // 'sub-region'
            });

            origins[0] = mesh_.cellCentres()[celli];
            distSqrs[0] = 0.0;  // Initial distance
            surfaces[0] = zoneTag;
            allCellInfo[celli] = wallPoints(origins, distSqrs, surfaces);
        }
    }


    DynamicList<wallPoints> faceDist(mesh_.nFaces()/128);
    DynamicList<label> changedFaces(mesh_.nFaces()/128);

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const label own = mesh_.faceOwner()[facei];
        const label nei = mesh_.faceNeighbour()[facei];
        if (cellToZone[own] >= 0 && cellToZone[nei] < 0)
        {
            // boundary between cellZone (own) and background/unvisited (nei)

            origins[0] = mesh_.faceCentres()[facei];
            distSqrs[0] = 0.0;  // Initial distance
            surfaces[0] = FixedList<label, 3>
            ({
                cellToZone[own],        // zone
                0,                      // 'region'
                labelMax                // 'sub-region'
            });
            faceDist.append(wallPoints(origins, distSqrs, surfaces));
            changedFaces.append(facei);
        }
        else if (cellToZone[own] < 0 && cellToZone[nei] >= 0)
        {
            // boundary between cellZone (nei) and background/unvisited (own)

            origins[0] = mesh_.faceCentres()[facei];
            distSqrs[0] = 0.0;  // Initial distance
            surfaces[0] = FixedList<label, 3>
            ({
                cellToZone[nei],        // zone
                0,                      // 'region'
                labelMax                // 'sub-region'
            });
            faceDist.append(wallPoints(origins, distSqrs, surfaces));
            changedFaces.append(facei);
        }
        else if
        (
            unnamedSurfaceRegion1[facei] != -1
         || unnamedSurfaceRegion2[facei] != -1
        )
        {
            // Seed (yet unpatched) wall faces

            origins[0] = mesh_.faceCentres()[facei];
            distSqrs[0] = 0.0;  // Initial distance
            surfaces[0] = wallTag;
            faceDist.append(wallPoints(origins, distSqrs, surfaces));
            changedFaces.append(facei);
        }
    }

    // Get coupled neighbour cellRegion
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Calculate region to zone from cellRegions on either side of coupled
    // face.
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            // Like internal faces
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label own = mesh_.faceOwner()[facei];
                label bFacei = facei-mesh_.nInternalFaces();
                if (cellToZone[own] >= 0 && neiCellZone[bFacei] < 0)
                {
                    origins[0] = mesh_.faceCentres()[facei];
                    distSqrs[0] = 0.0;
                    surfaces[0] = FixedList<label, 3>
                    ({
                        cellToZone[own],        // zone
                        0,                      // 'region'
                        labelMax                // 'sub-region'
                    });
                    faceDist.append(wallPoints(origins, distSqrs, surfaces));
                    changedFaces.append(facei);
                }
                else if (cellToZone[own] < 0 && neiCellZone[bFacei] >= 0)
                {
                    // Handled on nbr processor
                }
                else if
                (
                    unnamedSurfaceRegion1[facei] != -1
                 || unnamedSurfaceRegion2[facei] != -1
                )
                {
                    // Seed (yet unpatched) wall faces
                    origins[0] = mesh_.faceCentres()[facei];
                    distSqrs[0] = 0.0;  // Initial distance
                    surfaces[0] = wallTag;
                    faceDist.append(wallPoints(origins, distSqrs, surfaces));
                    changedFaces.append(facei);
                }
            }
        }
        else
        {
            // Seed wall faces
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label own = mesh_.faceOwner()[facei];
                if (cellToZone[own] < 0)
                {
                    origins[0] = mesh_.faceCentres()[facei];
                    distSqrs[0] = 0.0;  // Initial distance
                    surfaces[0] = wallTag;
                    faceDist.append(wallPoints(origins, distSqrs, surfaces));
                    changedFaces.append(facei);
                }
            }
        }
    }


    List<wallPoints> allFaceInfo(mesh_.nFaces());

    // No blocked faces, limitless gap size
    const bitSet isBlockedFace(mesh_.nFaces());
    //List<scalarList> regionToBlockSize(surfaces_.surfaces().size());
    //{
    //    forAll(surfaces_.surfaces(), surfi)
    //    {
    //        const label geomi = surfaces_.surfaces()[surfi];
    //        const auto& s = surfaces_.geometry()[geomi];
    //        const label nRegions = s.regions().size();
    //        regionToBlockSize[surfi].setSize(nRegions, Foam::sqr(GREAT));
    //    }
    //}

    //- regionToBlockSize is indexed with cellZone index (>= 0) and region
    //- on cellZone (currently always 0)
    const auto& czs = mesh_.cellZones();
    List<scalarList> regionToBlockSize(czs.size());
    for (auto& blockSizes : regionToBlockSize)
    {
        blockSizes.setSize(1, Foam::sqr(GREAT));
    }


    wallPoints::trackData td(isBlockedFace, regionToBlockSize);
    FaceCellWave<wallPoints, wallPoints::trackData> wallDistCalc
    (
        mesh_,
        changedFaces,
        faceDist,
        allFaceInfo,
        allCellInfo,
        0,            // max iterations
        td
    );
    wallDistCalc.iterate(nGrowCellZones);


    // Check for cells which are within nGrowCellZones of two cellZones or
    // one cellZone and a wall
    // TBD. Currently only one layer of growth handled.

    bitSet isChangedCell(mesh_.nCells());

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].valid(wallDistCalc.data()))
        {
            const List<FixedList<label, 3>>& surfaces =
                allCellInfo[celli].surface();

            if (surfaces.size())
            {
                // Cell close to cellZone. Remove any free-standing baffles.
                // Done by marking as changed cell. Wip.
                isChangedCell.set(celli);
            }

            if (surfaces.size() > 1)
            {
                // Check if inbetween two cellZones or cellZone and wall
                // Find 'nearest' other cellZone
                scalar minDistSqr = Foam::sqr(GREAT);
                label minZone = -1;
                for (label i = 0; i < surfaces.size(); i++)
                {
                    const label zonei = surfaces[i][0];
                    const scalar d2 = allCellInfo[celli].distSqr()[i];
                    if (zonei != cellToZone[celli] && d2 < minDistSqr)
                    {
                        minDistSqr = d2;
                        minZone = zonei;   // zoneID
                    }
                }

                if (minDistSqr < Foam::sqr(GREAT))
                {
                    if (minZone != cellToZone[celli] && minZone != wallTag[0])
                    {
                        cellToZone[celli] = minZone;
                        isChangedCell.set(celli);
                    }
                }
            }
        }
    }

    // Fix any left-over unvisited cells
    if (backgroundZoneID != -2)
    {
        forAll(cellToZone, celli)
        {
            if (cellToZone[celli] == -2)
            {
                cellToZone[celli] = backgroundZoneID;
            }
        }
    }


    // Make sure to unset faces of changed cell

    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);


    label nUnnamed = 0;
    label nNamed = 0;
    for (const label celli : isChangedCell)
    {
        const cell& cFaces = mesh_.cells()[celli];
        for (const label facei : cFaces)
        {
            const label own = mesh_.faceOwner()[facei];
            const label ownZone = cellToZone[own];

            label nbrZone = -1;
            if (mesh_.isInternalFace(facei))
            {
                const label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
                nbrZone = (own != celli ? ownZone : neiZone);
            }
            else
            {
                nbrZone = neiCellZone[facei-mesh_.nInternalFaces()];
            }

            if (nbrZone == cellToZone[celli])
            {
                if
                (
                    unnamedSurfaceRegion1[facei] != -1
                 || unnamedSurfaceRegion2[facei] != -1
                )
                {
                    unnamedSurfaceRegion1[facei] = -1;
                    unnamedSurfaceRegion2[facei] = -1;
                    nUnnamed++;
                }
                if
                (
                    namedSurfaceRegion.size()
                 && namedSurfaceRegion[facei] != -1
                )
                {
                    namedSurfaceRegion[facei] = -1;
                    nNamed++;
                }
            }
        }
    }

    reduce(nUnnamed, sumOp<label>());
    reduce(nNamed, sumOp<label>());

    // Do always; might bypass if nNamed,nUnnamed zero
    syncTools::syncFaceList
    (
        mesh_,
        unnamedSurfaceRegion1,
        maxEqOp<label>()
    );
    syncTools::syncFaceList
    (
        mesh_,
        unnamedSurfaceRegion2,
        maxEqOp<label>()
    );
    if (namedSurfaceRegion.size())
    {
        syncTools::syncFaceList
        (
            mesh_,
            namedSurfaceRegion,
            maxEqOp<label>()
        );
    }

    if (debug)
    {
        Pout<< "growCellZone : grown cellZones by "
            << returnReduce(isChangedCell.count(), sumOp<label>())
            << " cells (moved from background to nearest cellZone)"
            << endl;
        Pout<< "growCellZone : unmarked " << nUnnamed
            << " unzoned intersections; " << nNamed << " zoned intersections; "
            << endl;
    }
}


void Foam::meshRefinement::makeConsistentFaceIndex
(
    const labelList& surfaceMap,
    const labelList& cellToZone,
    labelList& namedSurfaceRegion
) const
{
    // Make namedSurfaceRegion consistent with cellToZone - clear out any
    // blocked faces inbetween same cell zone (or background (=-1))
    // Do not do this for surfaces relating to 'pure' faceZones i.e.
    // faceZones without a cellZone. Note that we cannot check here
    // for different cellZones on either side but no namedSurfaceRegion
    // since cellZones can now originate from locationsInMesh as well
    // (instead of only through named surfaces)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label ownZone = cellToZone[faceOwner[faceI]];
        label neiZone = cellToZone[faceNeighbour[faceI]];
        label globalI = namedSurfaceRegion[faceI];

        if
        (
            ownZone == neiZone
         && globalI != -1
         && surfaceMap[surfaces_.whichSurface(globalI).first()] == -1
        )
        {
            namedSurfaceRegion[faceI] = -1;
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get coupled neighbour cellZone
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    // Use coupled cellZone to do check
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                label ownZone = cellToZone[faceOwner[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];
                label globalI = namedSurfaceRegion[faceI];

                if
                (
                    ownZone == neiZone
                 && globalI != -1
                 && surfaceMap[surfaces_.whichSurface(globalI).first()] == -1
                )
                {
                    namedSurfaceRegion[faceI] = -1;
                }
            }
        }
        else
        {
            // Unzonify boundary faces
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                label globalI = namedSurfaceRegion[faceI];

                if
                (
                    globalI != -1
                 && surfaceMap[surfaces_.whichSurface(globalI).first()] == -1
                )
                {
                    namedSurfaceRegion[faceI] = -1;
                }
            }
        }
    }
}


void Foam::meshRefinement::getIntersections
(
    const labelList& surfacesToTest,
    const pointField& neiCc,
    const labelList& testFaces,

    labelList& namedSurfaceRegion,
    bitSet& posOrientation
) const
{
    namedSurfaceRegion.setSize(mesh_.nFaces());
    namedSurfaceRegion = -1;

    posOrientation.setSize(mesh_.nFaces());
    posOrientation = false;

    // Statistics: number of faces per faceZone
    labelList nSurfFaces(surfaces_.surfZones().size(), Zero);

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    {
        labelList minLevel;
        calcCellCellRays
        (
            neiCc,
            labelList(neiCc.size(), -1),
            testFaces,
            start,
            end,
            minLevel
        );
    }

    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note that we intersect all intersected faces again. Could reuse
    // the information already in surfaceIndex_.

    labelList surface1;
    labelList region1;
    List<pointIndexHit> hit1;
    vectorField normal1;
    labelList surface2;
    labelList region2;
    List<pointIndexHit> hit2;
    vectorField normal2;

    surfaces_.findNearestIntersection
    (
        surfacesToTest,
        start,
        end,

        surface1,
        hit1,
        region1,
        normal1,

        surface2,
        hit2,
        region2,
        normal2
    );

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];
        const vector& area = mesh_.faceAreas()[faceI];

        if (surface1[i] != -1)
        {
            // If both hit should probably choose 'nearest'
            if
            (
                surface2[i] != -1
             && (
                    magSqr(hit2[i].hitPoint())
                  < magSqr(hit1[i].hitPoint())
                )
            )
            {
                namedSurfaceRegion[faceI] = surfaces_.globalRegion
                (
                    surface2[i],
                    region2[i]
                );
                posOrientation.set(faceI, ((area&normal2[i]) > 0));
                nSurfFaces[surface2[i]]++;
            }
            else
            {
                namedSurfaceRegion[faceI] = surfaces_.globalRegion
                (
                    surface1[i],
                    region1[i]
                );
                posOrientation.set(faceI, ((area&normal1[i]) > 0));
                nSurfFaces[surface1[i]]++;
            }
        }
        else if (surface2[i] != -1)
        {
            namedSurfaceRegion[faceI] = surfaces_.globalRegion
            (
                surface2[i],
                region2[i]
            );
            posOrientation.set(faceI, ((area&normal2[i]) > 0));
            nSurfFaces[surface2[i]]++;
        }
    }


    // surfaceIndex might have different surfaces on both sides if
    // there happen to be a (obviously thin) surface with different
    // regions between the cell centres. If one is on a named surface
    // and the other is not this might give problems so sync.
    syncTools::syncFaceList
    (
        mesh_,
        namedSurfaceRegion,
        maxEqOp<label>()
    );

    // Print a bit
    if (debug)
    {
        forAll(nSurfFaces, surfI)
        {
            Pout<< "Surface:"
                << surfaces_.names()[surfI]
                << "  nZoneFaces:" << nSurfFaces[surfI] << nl;
        }
        Pout<< endl;
    }
}


void Foam::meshRefinement::zonify
(
    const bool allowFreeStandingZoneFaces,
    const label nErodeCellZones,
    const label backgroundZoneID,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh,
    const bool exitIfLeakPath,
    const refPtr<coordSetWriter>& leakPathFormatter,
    refPtr<surfaceWriter>& surfFormatter,

    labelList& cellToZone,
    labelList& unnamedRegion1,
    labelList& unnamedRegion2,
    labelList& namedSurfaceRegion,
    bitSet& posOrientation
) const
{
    // Determine zones for cells and faces
    // cellToZone:
    // -2  : unset
    // -1  : not in any zone (zone 'none' or background zone)
    // >=0 : zoneID
    // namedSurfaceRegion, posOrientation:
    // -1  : face not intersected by named surface
    // >=0 : globalRegion (surface+region)
    //       (and posOrientation: surface normal v.s. face normal)

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();

    // Collect inside and outside into single list
    const List<pointField> allLocations
    (
        refinementParameters::zonePoints
        (
            locationsInMesh,
            zonesInMesh,
            locationsOutsideMesh
        )
    );

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);


    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));
    labelList unnamedSurfaces(surfaceZonesInfo::getUnnamedSurfaces(surfZones));

    // Get map from surface to cellZone (or -1)
    labelList surfaceToCellZone;
    if (namedSurfaces.size())
    {
        // Get/add cellZones corresponding to surface names
        surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
    }


    cellToZone.setSize(mesh_.nCells());
    cellToZone = -2;

    namedSurfaceRegion.clear();
    posOrientation.clear();



    // 1. Test all (unnamed & named) surfaces

    // Unnamed surfaces. Global surface region of intersection (or -1)
    getIntersections
    (
        unnamedSurfaces,
        neiCc,
        intersectedFaces(),
        unnamedRegion1,
        unnamedRegion2
    );

    // Extend with hole closing faces (only if locationsOutsideMesh)
    labelList unnamedFaces;
    labelList unnamedClosureFaces;
    labelList unnamedToClosure;
    autoPtr<mapDistribute> unnamedMapPtr;
    if (locationsOutsideMesh.size())
    {
        unnamedFaces = ListOps::findIndices
        (
            unnamedRegion1,
            [](const label x){return x != -1;}
        );

        const globalIndex globalUnnamedFaces(unnamedFaces.size());

        unnamedMapPtr = holeToFace::calcClosure
        (
            mesh_,
            allLocations,
            unnamedFaces,
            globalUnnamedFaces,
            true,                   // allow erosion

            unnamedClosureFaces,
            unnamedToClosure
        );

        if (debug)
        {
            Pout<< "meshRefinement::zonify : found wall closure faces:"
                << unnamedClosureFaces.size()
                << "  map:" << bool(unnamedMapPtr) << endl;
        }


        // Add to unnamedRegion1, unnamedRegion2
        if (unnamedMapPtr)
        {
            // Dump leak path
            if (leakPathFormatter)
            {
                boolList blockedFace(mesh_.nFaces(), false);
                UIndirectList<bool>(blockedFace, unnamedFaces) = true;
                const fileName fName
                (
                    writeLeakPath
                    (
                        mesh_,
                        locationsInMesh,
                        locationsOutsideMesh,
                        blockedFace,
                        leakPathFormatter.constCast()
                    )
                );
                Info<< "Dumped leak path to " << fName << endl;
            }

            auto& err =
            (
                exitIfLeakPath
              ? FatalErrorInFunction
              : WarningInFunction
            );

            err << "Locations in mesh " << locationsInMesh
                << " connect to one of the locations outside mesh "
                << locationsOutsideMesh << endl;

            if (exitIfLeakPath)
            {
                FatalError << exit(FatalError);
            }


            labelList packedRegion1
            (
                UIndirectList<label>(unnamedRegion1, unnamedFaces)
            );
            unnamedMapPtr->distribute(packedRegion1);
            labelList packedRegion2
            (
                UIndirectList<label>(unnamedRegion2, unnamedFaces)
            );
            unnamedMapPtr->distribute(packedRegion2);
            forAll(unnamedClosureFaces, i)
            {
                const label sloti = unnamedToClosure[i];

                if (sloti != -1)
                {
                    const label facei = unnamedClosureFaces[i];
                    const label region1 = unnamedRegion1[facei];
                    const label slotRegion1 = packedRegion1[sloti];
                    const label region2 = unnamedRegion2[facei];
                    const label slotRegion2 = packedRegion2[sloti];

                    if (slotRegion1 != region1 || slotRegion2 != region2)
                    {
                        unnamedRegion1[facei] = slotRegion1;
                        unnamedRegion2[facei] = slotRegion2;
                    }
                }
            }
        }
    }


    // Extend with hole closing faces (only if locationsOutsideMesh)
    labelList namedFaces;
    labelList namedClosureFaces;
    labelList namedToClosure;
    autoPtr<mapDistribute> namedMapPtr;
    if (namedSurfaces.size())
    {
        getIntersections
        (
            namedSurfaces,
            neiCc,
            intersectedFaces(),
            namedSurfaceRegion,
            posOrientation
        );

        // Ideally we'd like to close 'cellZone' surfaces. The problem is
        // that we don't (easily) know which locationsInMesh should be inside
        // the surface and which aren't. With 'insidePoint' definition of
        // cellZone we have a location inside the cellZone but how do we
        // know where the locationsInMesh are? Are they inside the cellZone
        // as well? Only with the 'locationsInMesh' notation where we specify
        // the cellZone and the seedpoint could we make sure that we cannot
        // walk from one to the other.
        // For now disable hole closure on cellZones

        //if (locationsOutsideMesh.size())
        //{
        //    namedFaces = ListOps::findIndices
        //    (
        //        namedSurfaceRegion,
        //        [](const label x){return x != -1;}
        //    );
        //
        //    {
        //        OBJstream str(mesh_.time().timePath()/"namedFaces.obj");
        //        Pout<< "Writing " << namedFaces.size() << " zone faces to "
        //            << str.name() << endl;
        //        str.write
        //        (
        //            UIndirectList<face>(mesh_.faces(), namedFaces)(),
        //            mesh_.points()
        //        );
        //    }
        //
        //
        //    const globalIndex globalNamedFaces(namedFaces.size());
        //
        //    namedMapPtr = holeToFace::calcClosure
        //    (
        //        mesh_,
        //        allLocations,
        //        namedFaces,             // or also unnamedFaces?
        //        globalNamedFaces,
        //        true,                   // allow erosion
        //
        //        namedClosureFaces,
        //        namedToClosure
        //    );
        //
        //    if (debug)
        //    {
        //        Pout<< "meshRefinement::zonify :"
        //            << " found faceZone closure faces:"
        //            << namedClosureFaces.size()
        //            << "  map:" << bool(namedMapPtr) << endl;
        //    }
        //
        //    // Add to namedSurfaceRegion, posOrientation
        //    if (namedMapPtr)
        //    {
        //        WarningInFunction
        //            << "Detected and closed leak path"
        //            << " through zoned surfaces from "
        //            << locationsInMesh << " to " << locationsOutsideMesh
        //            << endl;
        //
        //        // Dump leak path
        //        if (leakPathFormatter)
        //        {
        //            boolList blockedFace(mesh_.nFaces(), false);
        //            UIndirectList<bool>(blockedFace, unnamedFaces) = true;
        //            UIndirectList<bool>(blockedFace, namedFaces) = true;
        //            const fileName fName
        //            (
        //                writeLeakPath
        //                (
        //                    mesh_,
        //                    locationsInMesh,
        //                    locationsOutsideMesh,
        //                    blockedFace,
        //                    leakPathFormatter.constCast()
        //                )
        //            );
        //            Info<< "Dumped leak path to " << fName << endl;
        //        }
        //
        //        labelList packedSurfaceRegion
        //        (
        //            UIndirectList<label>(namedSurfaceRegion, namedFaces)
        //        );
        //        namedMapPtr->distribute(packedSurfaceRegion);
        //        boolList packedOrientation(posOrientation.size());
        //        forAll(namedFaces, i)
        //        {
        //            const label facei = namedFaces[i];
        //            packedOrientation[i] = posOrientation[facei];
        //        }
        //        namedMapPtr->distribute(packedOrientation);
        //        forAll(namedClosureFaces, i)
        //        {
        //            const label sloti = namedToClosure[i];
        //            if (sloti != -1)
        //            {
        //                const label facei = namedClosureFaces[i];
        //                const label regioni = namedSurfaceRegion[facei];
        //                const label slotRegioni = packedSurfaceRegion[sloti];
        //                const bool orient = posOrientation[facei];
        //                const bool slotOrient = packedOrientation[sloti];
        //
        //                if (slotRegioni != regioni || slotOrient != orient)
        //                {
        //                    namedSurfaceRegion[facei] = slotRegioni;
        //                    posOrientation[facei] = slotOrient;
        //                }
        //            }
        //        }
        //    }
        //}
    }


    // 1b. Add any hole closure faces to frozenPoints pointZone
    {
        bitSet isClosureFace(mesh_.nFaces());
        isClosureFace.set(unnamedClosureFaces);
        isClosureFace.set(namedClosureFaces);
        const labelList closureFaces(isClosureFace.sortedToc());

        const uindirectPrimitivePatch pp
        (
            UIndirectList<face>(mesh_.faces(), closureFaces),
            mesh_.points()
        );

        // Count number of faces per edge
        const labelList nEdgeFaces(countEdgeFaces(pp));

        // Freeze all internal points
        bitSet isFrozenPoint(mesh_.nPoints());
        forAll(nEdgeFaces, edgei)
        {
            if (nEdgeFaces[edgei] != 1)
            {
                const edge& e = pp.edges()[edgei];
                isFrozenPoint.set(pp.meshPoints()[e[0]]);
                isFrozenPoint.set(pp.meshPoints()[e[1]]);
            }
        }

        // Lookup/add pointZone and include its points
        pointZoneMesh& pointZones =
            const_cast<pointZoneMesh&>(mesh_.pointZones());
        const label zonei =
            const_cast<meshRefinement&>(*this).addPointZone("frozenPoints");
        const bitSet oldSet(mesh_.nPoints(), pointZones[zonei]);
        isFrozenPoint.set(oldSet);

        syncTools::syncPointList
        (
            mesh_,
            isFrozenPoint,
            orEqOp<unsigned int>(),
            0u
        );

        // Override addressing
        pointZones.clearAddressing();
        pointZones[zonei] = isFrozenPoint.sortedToc();

        if (debug)
        {
            mkDir(mesh_.time().timePath());
            const pointZone& pz = pointZones[zonei];
            OBJstream str(mesh_.time().timePath()/pz.name()+".obj");
            Pout<< "Writing " << pz.size() << " frozen points to "
                << str.name() << endl;
            for (const label pointi : pz)
            {
                str.write(mesh_.points()[pointi]);
            }
        }

        if (returnReduceOr(unnamedClosureFaces.size()) && surfFormatter)
        {
            fileName outputName
            (
                mesh_.time().globalPath()
              / functionObject::outputPrefix
              / mesh_.pointsInstance()
              / "unnamedClosureFaces"
            );
            outputName.clean();  // Remove unneeded ".."

            Info<< "Writing "
                << returnReduce(unnamedClosureFaces.size(), sumOp<label>())
                << " unnamedClosureFaces to " << outputName << endl;

            const indirectPrimitivePatch setPatch
            (
                IndirectList<face>(mesh_.faces(), unnamedClosureFaces),
                mesh_.points()
            );

            surfFormatter->open
            (
                setPatch.localPoints(),
                setPatch.localFaces(),
                outputName
            );
            surfFormatter->write();
            surfFormatter->clear();
        }
        if (returnReduceOr(namedClosureFaces.size()) && surfFormatter)
        {
            fileName outputName
            (
                mesh_.time().globalPath()
              / functionObject::outputPrefix
              / mesh_.pointsInstance()
              / "namedClosureFaces"
            );
            outputName.clean();  // Remove unneeded ".."

            Info<< "Writing "
                << returnReduce(namedClosureFaces.size(), sumOp<label>())
                << " namedClosureFaces to " << outputName << endl;

            const indirectPrimitivePatch setPatch
            (
                IndirectList<face>(mesh_.faces(), namedClosureFaces),
                mesh_.points()
            );

            surfFormatter->open
            (
                setPatch.localPoints(),
                setPatch.localFaces(),
                outputName
            );
            surfFormatter->write();
            surfFormatter->clear();
        }
    }



    // 2. Walk from locationsInMesh. Hard set cellZones.
    //    Note: walk through faceZones! (these might get overridden later)

    if (locationsInMesh.size())
    {
        Info<< "Setting cellZones according to locationsInMesh:" << endl;

        labelList locationsZoneIDs(zonesInMesh.size(), -1);
        forAll(locationsInMesh, i)
        {
            const word& name = zonesInMesh[i];

            Info<< "Location : " << locationsInMesh[i] << nl
                << "    cellZone : " << name << endl;

            if (name != "none")
            {
                label zoneID = mesh_.cellZones().findZoneID(name);
                if (zoneID == -1)
                {
                    FatalErrorInFunction << "problem" << abort(FatalError);
                }
                locationsZoneIDs[i] = zoneID;
            }
        }
        Info<< endl;


        // Assign cellZone according to seed points. Walk through named surfaces
        findCellZoneInsideWalk
        (
            locationsInMesh,    // locations
            locationsZoneIDs,   // index of cellZone (or -1)
            unnamedRegion1,     // per face -1 (unblocked) or >= 0 (blocked)
            cellToZone
        );
    }


    // 3. Mark named-surfaces-with-insidePoint. Hard set cellZones.

    labelList locationSurfaces
    (
        surfaceZonesInfo::getInsidePointNamedSurfaces(surfZones)
    );

    if (locationSurfaces.size())
    {
        Info<< "Found " << locationSurfaces.size()
            << " named surfaces with a provided inside point."
            << " Assigning cells inside these surfaces"
            << " to the corresponding cellZone."
            << nl << endl;

        // Collect per surface the -insidePoint -cellZone
        // Usually only a single inside point per surface so no clever
        // counting - just use DynamicField
        DynamicField<point> insidePoints(locationSurfaces.size());
        DynamicList<label> insidePointCellZoneIDs(locationSurfaces.size());
        forAll(locationSurfaces, i)
        {
            const label surfI = locationSurfaces[i];
            const auto& surfInsidePoints = surfZones[surfI].zoneInsidePoints();

            const word& name = surfZones[surfI].cellZoneName();
            label zoneID = -1;
            if (name != "none")
            {
                zoneID = mesh_.cellZones().findZoneID(name);
                if (zoneID == -1)
                {
                    FatalErrorInFunction
                        << "Specified non-existing cellZone " << name
                        << " for surface " << surfaces_.names()[surfI]
                        << abort(FatalError);
                }
            }

            for (const auto& pt : surfInsidePoints)
            {
                insidePoints.append(pt);
                insidePointCellZoneIDs.append(zoneID);
            }
        }


        // Stop at unnamed or named surface
        labelList allRegion1(mesh_.nFaces(), -1);
        forAll(namedSurfaceRegion, faceI)
        {
            allRegion1[faceI] = max
            (
                unnamedRegion1[faceI],
                namedSurfaceRegion[faceI]
            );
        }

        findCellZoneInsideWalk
        (
            insidePoints,           // locations
            insidePointCellZoneIDs, // index of cellZone
            allRegion1,             // per face -1 (unblocked) or >= 0 (blocked)
            cellToZone
        );
    }


    // 4. Mark named-surfaces-with-geometric faces. Do geometric test. Soft set
    // cellZones. Correct through making consistent.

    // Closed surfaces with cellZone specified.
    labelList closedNamedSurfaces
    (
        surfaceZonesInfo::getClosedNamedSurfaces
        (
            surfZones,
            surfaces_.geometry(),
            surfaces_.surfaces()
        )
    );

    if (closedNamedSurfaces.size())
    {
        Info<< "Found " << closedNamedSurfaces.size()
            << " closed, named surfaces. Assigning cells in/outside"
            << " these surfaces to the corresponding cellZone."
            << nl << endl;

        findCellZoneGeometric
        (
            neiCc,
            closedNamedSurfaces,    // indices of closed surfaces
            namedSurfaceRegion,     // per face index of named surface + region
            surfaceToCellZone,      // cell zone index per surface

            cellToZone
        );
    }


    // 5. Find any unassigned regions (from regionSplit)

    if (namedSurfaces.size())
    {
        if (nErodeCellZones == 0)
        {
            Info<< "Walking from known cellZones; crossing a faceZone "
                << "face changes cellZone" << nl << endl;

            // Put unassigned regions into any connected cellZone
            findCellZoneTopo
            (
                backgroundZoneID,
                pointField(0),
                unnamedRegion1,      // Intersections with unnamed surfaces
                namedSurfaceRegion,  // Intersections with named surfaces
                surfaceToCellZone,
                cellToZone
            );
        }
        else if (nErodeCellZones < 0)
        {
            Info<< "Growing cellZones by " << -nErodeCellZones
                << " layers" << nl << endl;

            growCellZone
            (
                -nErodeCellZones,
                backgroundZoneID,
                unnamedRegion1,
                unnamedRegion2,
                namedSurfaceRegion,
                cellToZone
            );
        }
        else
        {
            Info<< "Eroding cellZone cells to make these consistent with"
                << " faceZone faces" << nl << endl;

            // Erode cell zone cells (connected to an unassigned cell)
            // and put them into backgroundZone
            erodeCellZone
            (
                nErodeCellZones,
                backgroundZoneID,
                unnamedRegion1,
                namedSurfaceRegion,
                cellToZone
            );
        }


        // Make sure namedSurfaceRegion is unset inbetween same cell zones.
        if (!allowFreeStandingZoneFaces)
        {
            Info<< "Only keeping zone faces inbetween different cellZones."
                << nl << endl;

            // Surfaces with faceZone but no cellZone
            labelList standaloneNamedSurfaces
            (
                surfaceZonesInfo::getStandaloneNamedSurfaces
                (
                    surfZones
                )
            );

            // Construct map from surface index to index in
            // standaloneNamedSurfaces (or -1)
            labelList surfaceMap(surfZones.size(), -1);
            forAll(standaloneNamedSurfaces, i)
            {
                surfaceMap[standaloneNamedSurfaces[i]] = i;
            }

            makeConsistentFaceIndex
            (
                surfaceMap,
                cellToZone,
                namedSurfaceRegion
            );
        }
    }
    else if (nErodeCellZones < 0 && gMax(cellToZone) >= 0)
    {
        // With multiple locationsInMesh and non-trivial cellZones we allow
        // some growing of the cellZones to avoid any background cells

        Info<< "Growing cellZones by " << -nErodeCellZones
            << " layers" << nl << endl;

        growCellZone
        (
            -nErodeCellZones,
            backgroundZoneID,
            unnamedRegion1,
            unnamedRegion2,
            namedSurfaceRegion, // note: potentially zero sized
            cellToZone
        );

        // Make sure namedSurfaceRegion is unset inbetween same cell zones.
        if (!allowFreeStandingZoneFaces && namedSurfaceRegion.size())
        {
            Info<< "Only keeping zone faces inbetween different cellZones."
                << nl << endl;

            // Surfaces with faceZone but no cellZone
            labelList standaloneNamedSurfaces
            (
                surfaceZonesInfo::getStandaloneNamedSurfaces
                (
                    surfZones
                )
            );

            // Construct map from surface index to index in
            // standaloneNamedSurfaces (or -1)
            labelList surfaceMap(surfZones.size(), -1);
            forAll(standaloneNamedSurfaces, i)
            {
                surfaceMap[standaloneNamedSurfaces[i]] = i;
            }

            makeConsistentFaceIndex
            (
                surfaceMap,
                cellToZone,
                namedSurfaceRegion
            );
        }
    }


    // Some stats
    if (debug)
    {
        label nZones = gMax(cellToZone)+1;

        label nUnvisited = 0;
        label nBackgroundCells = 0;
        labelList nZoneCells(nZones, Zero);
        forAll(cellToZone, cellI)
        {
            label zoneI = cellToZone[cellI];
            if (zoneI >= 0)
            {
                nZoneCells[zoneI]++;
            }
            else if (zoneI == -1)
            {
                nBackgroundCells++;
            }
            else if (zoneI == -2)
            {
                nUnvisited++;
            }
            else
            {
                FatalErrorInFunction
                    << "problem" << exit(FatalError);
            }
        }
        reduce(nUnvisited, sumOp<label>());
        reduce(nBackgroundCells, sumOp<label>());
        forAll(nZoneCells, zoneI)
        {
            reduce(nZoneCells[zoneI], sumOp<label>());
        }
        Info<< "nUnvisited      :" << nUnvisited << endl;
        Info<< "nBackgroundCells:" << nBackgroundCells << endl;
        Info<< "nZoneCells      :" << nZoneCells << endl;
    }
    if (debug&MESH)
    {
        const_cast<Time&>(mesh_.time())++;
        Pout<< "Writing cell zone allocation on mesh to time "
            << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            mesh_.time().path()/"cell2Zone"
        );
        volScalarField volCellToZone
        (
            IOobject
            (
                "cellToZone",
                timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        forAll(cellToZone, cellI)
        {
            volCellToZone[cellI] = cellToZone[cellI];
        }
        volCellToZone.write();


        //mkDir(mesh_.time().path()/timeName());
        //OBJstream str
        //(
        //    mesh_.time().path()/timeName()/"zoneBoundaryFaces.obj"
        //);
        //Pout<< "Writing zone boundaries to " << str.name() << endl;
        //for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        //{
        //    const label ownZone = cellToZone[mesh_.faceOwner()[facei]];
        //    const label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
        //    if (ownZone != neiZone)
        //    {
        //        str.write(mesh_.faces()[facei], mesh_.points(), false);
        //    }
        //}
        //labelList neiCellZone;
        //syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
        //for
        //(
        //    label facei = mesh_.nInternalFaces();
        //    facei < mesh_.nFaces();
        //    facei++
        //)
        //{
        //    const label ownZone = cellToZone[mesh_.faceOwner()[facei]];
        //    const label bFacei = facei-mesh_.nInternalFaces();
        //    const label neiZone = neiCellZone[bFacei];
        //    if (ownZone != neiZone)
        //    {
        //        str.write(mesh_.faces()[facei], mesh_.points(), false);
        //    }
        //}

        //mkDir(mesh_.time().path()/timeName());
        //OBJstream str1
        //(
        //    mesh_.time().path()/timeName()/"unnamedRegion1.obj"
        //);
        //OBJstream str2
        //(
        //    mesh_.time().path()/timeName()/"unnamedRegion2.obj"
        //);
        //Pout<< "Writing unnamed1 to " << str1.name() << endl;
        //Pout<< "Writing unnamed2 to " << str2.name() << endl;
        //for (label facei = 0; facei < mesh_.nFaces(); facei++)
        //{
        //    if
        //    (
        //        unnamedRegion1[facei] < -1
        //     || unnamedRegion2[facei] < -1
        //    )
        //    {
        //        FatalErrorInFunction << "face:"
        //            << mesh_.faceCentres()[facei]
        //            << " unnamed1:" << unnamedRegion1[facei]
        //            << " unnamed2:" << unnamedRegion2[facei]
        //            << exit(FatalError);
        //    }
        //
        //    if (unnamedRegion1[facei] >= 0)
        //    {
        //        str1.write(mesh_.faces()[facei], mesh_.points(), false);
        //    }
        //
        //    if (unnamedRegion2[facei] >= 0)
        //    {
        //        str2.write(mesh_.faces()[facei], mesh_.points(), false);
        //    }
        //}

        //if (namedSurfaceRegion.size())
        //{
        //    OBJstream strNamed
        //    (
        //        mesh_.time().path()/timeName()/"namedSurfaceRegion.obj"
        //    );
        //    Pout<< "Writing named to " << strNamed.name() << endl;
        //    for (label facei = 0; facei < mesh_.nFaces(); facei++)
        //    {
        //        const face& f = mesh_.faces()[facei];
        //        if (namedSurfaceRegion[facei] < -1)
        //        {
        //            FatalErrorInFunction << "face:"
        //                << mesh_.faceCentres()[facei]
        //                << " unnamed1:" << unnamedRegion1[facei]
        //                << " unnamed2:" << unnamedRegion2[facei]
        //                << " named:" << namedSurfaceRegion[facei]
        //                << exit(FatalError);
        //        }
        //        if (namedSurfaceRegion[facei] >= 0)
        //        {
        //            strNamed.write(f, mesh_.points(), false);
        //        }
        //    }
        //}
    }
}


void Foam::meshRefinement::handleSnapProblems
(
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
{
    Info<< nl
        << "Introducing baffles to block off problem cells" << nl
        << "----------------------------------------------" << nl
        << endl;

    labelList facePatch;
    labelList faceZone;
    if (useTopologicalSnapDetection)
    {
        markFacesOnProblemCells
        (
            motionDict,
            removeEdgeConnectedCells,
            perpendicularAngle,
            globalToMasterPatch,

            facePatch,
            faceZone
        );
    }
    else
    {
        markFacesOnProblemCellsGeometric
        (
            snapParams,
            motionDict,
            globalToMasterPatch,
            globalToSlavePatch,

            facePatch,
            faceZone
        );
    }
    Info<< "Analyzed problem cells in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug&MESH)
    {
        faceSet problemFaces(mesh_, "problemFaces", mesh_.nFaces()/100);

        forAll(facePatch, faceI)
        {
            if (facePatch[faceI] != -1)
            {
                problemFaces.insert(faceI);
            }
        }
        problemFaces.instance() = timeName();
        Pout<< "Dumping " << problemFaces.size()
            << " problem faces to " << problemFaces.objectPath() << endl;
        problemFaces.write();
    }

    Info<< "Introducing baffles to delete problem cells." << nl << endl;

    if (debug)
    {
        ++runTime;
    }


    // Add faces-to-baffle to faceZone. For now do this outside of topoChanges
    {
        const faceZoneMesh& fzs = mesh_.faceZones();
        List<DynamicList<label>> zoneToFaces(fzs.size());
        List<DynamicList<bool>> zoneToFlip(fzs.size());

        // Start off with original contents
        forAll(fzs, zonei)
        {
            zoneToFaces[zonei].append(fzs[zonei]);
            zoneToFlip[zonei].append(fzs[zonei].flipMap());
        }

        // Add any to-be-patched face
        forAll(facePatch, facei)
        {
            if (facePatch[facei] != -1)
            {
                label zonei = faceZone[facei];
                if (zonei != -1)
                {
                    zoneToFaces[zonei].append(facei);
                    zoneToFlip[zonei].append(false);
                }
            }
        }

        forAll(zoneToFaces, zonei)
        {
            surfaceZonesInfo::addFaceZone
            (
                fzs[zonei].name(),
                zoneToFaces[zonei],
                zoneToFlip[zonei],
                mesh_
            );
        }
    }


    // Create baffles with same owner and neighbour for now.
    createBaffles(facePatch, facePatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles", true);

    if (debug&MESH)
    {
        const_cast<Time&>(mesh_.time())++;
        Pout<< "Writing extra baffled mesh to time "
            << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"extraBaffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::labelList Foam::meshRefinement::freeStandingBaffleFaces
(
    const labelList& faceToZone,
    const labelList& cellToZone,
    const labelList& neiCellZone
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    // We want to pick up the faces to orient. These faces come in
    // two variants:
    // - faces originating from stand-alone faceZones
    //   (these will most likely have no cellZone on either side so
    //    ownZone and neiZone both -1)
    // - sticky-up faces originating from a 'bulge' in a outside of
    //   a cellZone. These will have the same cellZone on either side.
    //   How to orient these is not really clearly defined so do them
    //   same as stand-alone faceZone faces for now. (Normally these will
    //   already have been removed by the 'allowFreeStandingZoneFaces=false'
    //   default setting)

    // Note that argument neiCellZone will have -1 on uncoupled boundaries.

    DynamicList<label> faceLabels(mesh_.nFaces()/100);

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceToZone[faceI] != -1)
        {
            // Free standing baffle?
            label ownZone = cellToZone[faceOwner[faceI]];
            label neiZone = cellToZone[faceNeighbour[faceI]];
            if (ownZone == neiZone)
            {
                faceLabels.append(faceI);
            }
        }
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll(pp, i)
        {
            label faceI = pp.start()+i;
            if (faceToZone[faceI] != -1)
            {
                // Free standing baffle?
                label ownZone = cellToZone[faceOwner[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];
                if (ownZone == neiZone)
                {
                    faceLabels.append(faceI);
                }
            }
        }
    }
    return faceLabels.shrink();
}


void Foam::meshRefinement::calcPatchNumMasterFaces
(
    const bitSet& isMasterFace,
    const indirectPrimitivePatch& patch,
    labelList& nMasterFacesPerEdge
) const
{
    // Number of (master)faces per edge
    nMasterFacesPerEdge.setSize(patch.nEdges());
    nMasterFacesPerEdge = 0;

    forAll(patch.addressing(), faceI)
    {
        const label meshFaceI = patch.addressing()[faceI];

        if (isMasterFace.test(meshFaceI))
        {
            const labelList& fEdges = patch.faceEdges()[faceI];
            forAll(fEdges, fEdgeI)
            {
                nMasterFacesPerEdge[fEdges[fEdgeI]]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        patch.meshEdges(mesh_.edges(), mesh_.pointEdges()),
        nMasterFacesPerEdge,
        plusEqOp<label>(),
        label(0)
    );
}


Foam::label Foam::meshRefinement::markPatchZones
(
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    labelList& faceToZone
) const
{
    List<edgeTopoDistanceData<label>> allEdgeInfo(patch.nEdges());
    List<edgeTopoDistanceData<label>> allFaceInfo(patch.size());


    // Protect all non-manifold edges
    {
        // label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = edgeTopoDistanceData<label>(0, -2);
                // ++nProtected;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " non-manifold edges" << nl << endl;
    }


    // Hand out zones

    DynamicList<label> changedEdges;
    DynamicList<edgeTopoDistanceData<label>> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        edgeTopoDistanceData<label>
    >::propagationTol();

    int dummyTrackData;

    const globalIndex globalFaces(patch.size());

    label faceI = 0;

    label currentZoneI = 0;

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        for (; faceI < allFaceInfo.size(); faceI++)
        {
            if (!allFaceInfo[faceI].valid(dummyTrackData))
            {
                globalSeed = globalFaces.toGlobal(faceI);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        if (globalFaces.isLocal(globalSeed))
        {
            const label seedFaceI = globalFaces.toLocal(globalSeed);

            edgeTopoDistanceData<label>& faceInfo = allFaceInfo[seedFaceI];

            // Set face
            faceInfo = edgeTopoDistanceData<label>(0, currentZoneI);

            // .. and seed its edges
            const labelList& fEdges = patch.faceEdges()[seedFaceI];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                edgeTopoDistanceData<label>& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFaceI,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduceAnd(changedEdges.empty()))
        {
            break;
        }


        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            edgeTopoDistanceData<label>
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );

        currentZoneI++;
    }


    faceToZone.setSize(patch.size());
    forAll(allFaceInfo, faceI)
    {
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            FatalErrorInFunction
                << "Problem: unvisited face " << faceI
                << " at " << patch.faceCentres()[faceI]
                << exit(FatalError);
        }
        faceToZone[faceI] = allFaceInfo[faceI].data();
    }

    return currentZoneI;
}


void Foam::meshRefinement::consistentOrientation
(
    const bitSet& isMasterFace,
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    const labelList& faceToZone,
    const Map<label>& zoneToOrientation,
    bitSet& meshFlipMap
) const
{
    const polyBoundaryMesh& bm = mesh_.boundaryMesh();

    // Data on all edges and faces
    List<patchFaceOrientation> allEdgeInfo(patch.nEdges());
    List<patchFaceOrientation> allFaceInfo(patch.size());

    // Make sure we don't walk through
    // - slaves of coupled faces
    // - non-manifold edges
    {
        // label nProtected = 0;

        forAll(patch.addressing(), faceI)
        {
            const label meshFaceI = patch.addressing()[faceI];
            const label patchI = bm.whichPatch(meshFaceI);

            if
            (
                patchI != -1
             && bm[patchI].coupled()
             && !isMasterFace.test(meshFaceI)
            )
            {
                // Slave side. Mark so doesn't get visited.
                allFaceInfo[faceI] = orientedSurface::NOFLIP;
                // ++nProtected;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " slaves of coupled faces" << nl << endl;
    }
    {
        label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = orientedSurface::NOFLIP;
                ++nProtected;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " non-manifold edges" << nl << endl;
    }



    DynamicList<label> changedEdges;
    DynamicList<patchFaceOrientation> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchFaceOrientation
    >::propagationTol();

    int dummyTrackData;

    globalIndex globalFaces(patch.size());

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        forAll(allFaceInfo, faceI)
        {
            if (allFaceInfo[faceI] == orientedSurface::UNVISITED)
            {
                globalSeed = globalFaces.toGlobal(faceI);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        if (globalFaces.isLocal(globalSeed))
        {
            const label seedFaceI = globalFaces.toLocal(globalSeed);

            // Determine orientation of seedFace

            patchFaceOrientation& faceInfo = allFaceInfo[seedFaceI];

            // Start off with correct orientation
            faceInfo = orientedSurface::NOFLIP;

            if (zoneToOrientation[faceToZone[seedFaceI]] < 0)
            {
                faceInfo.flip();
            }


            const labelList& fEdges = patch.faceEdges()[seedFaceI];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchFaceOrientation& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFaceI,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduceAnd(changedEdges.empty()))
        {
            break;
        }



        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchFaceOrientation
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );
    }


    // Push master zone info over to slave (since slave faces never visited)
    {
        labelList neiStatus
        (
            mesh_.nBoundaryFaces(),
            orientedSurface::UNVISITED
        );

        forAll(patch.addressing(), i)
        {
            const label meshFaceI = patch.addressing()[i];
            if (!mesh_.isInternalFace(meshFaceI))
            {
                neiStatus[meshFaceI-mesh_.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiStatus);

        forAll(patch.addressing(), i)
        {
            const label meshFaceI = patch.addressing()[i];
            const label patchI = bm.whichPatch(meshFaceI);

            if
            (
                patchI != -1
             && bm[patchI].coupled()
             && !isMasterFace.test(meshFaceI)
            )
            {
                // Slave side. Take flipped from neighbour
                label bFaceI = meshFaceI-mesh_.nInternalFaces();

                if (neiStatus[bFaceI] == orientedSurface::NOFLIP)
                {
                    allFaceInfo[i] = orientedSurface::FLIP;
                }
                else if (neiStatus[bFaceI] == orientedSurface::FLIP)
                {
                    allFaceInfo[i] = orientedSurface::NOFLIP;
                }
                else
                {
                    FatalErrorInFunction
                        << "Incorrect status for face " << meshFaceI
                        << abort(FatalError);
                }
            }
        }
    }


    // Convert to meshFlipMap and adapt faceZones

    meshFlipMap.setSize(mesh_.nFaces());
    meshFlipMap = false;

    forAll(allFaceInfo, faceI)
    {
        label meshFaceI = patch.addressing()[faceI];

        if (allFaceInfo[faceI] == orientedSurface::NOFLIP)
        {
            meshFlipMap.unset(meshFaceI);
        }
        else if (allFaceInfo[faceI] == orientedSurface::FLIP)
        {
            meshFlipMap.set(meshFaceI);
        }
        else
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << faceI
                << " centre:" << mesh_.faceCentres()[meshFaceI]
                << abort(FatalError);
        }
    }
}


void Foam::meshRefinement::zonify
(
    // Get per face whether is it master (of a coupled set of faces)
    const bitSet& isMasterFace,
    const labelList& cellToZone,
    const labelList& neiCellZone,
    const labelList& faceToZone,
    const bitSet& meshFlipMap,
    polyTopoChange& meshMod
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label faceZoneI = faceToZone[faceI];

        if (faceZoneI != -1)
        {
            // Orient face zone to have slave cells in min cell zone.
            // Note: logic to use flipMap should be consistent with logic
            //       to pick up the freeStandingBaffleFaces!

            label ownZone = cellToZone[faceOwner[faceI]];
            label neiZone = cellToZone[faceNeighbour[faceI]];

            bool flip;

            if (ownZone == neiZone)
            {
                // free-standing face. Use geometrically derived orientation
                flip = meshFlipMap.test(faceI);
            }
            else
            {
                flip =
                (
                    ownZone == -1
                 || (neiZone != -1 && ownZone > neiZone)
                );
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[faceI],           // modified face
                    faceI,                          // label of face
                    faceOwner[faceI],               // owner
                    faceNeighbour[faceI],           // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    faceZoneI,                      // zone for face
                    flip                            // face flip in zone
                )
            );
        }
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Set owner as no-flip
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        label faceI = pp.start();

        forAll(pp, i)
        {
            label faceZoneI = faceToZone[faceI];

            if (faceZoneI != -1)
            {
                label ownZone = cellToZone[faceOwner[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];

                bool flip;

                if (ownZone == neiZone)
                {
                    // free-standing face. Use geometrically derived orientation
                    flip = meshFlipMap.test(faceI);
                }
                else
                {
                    flip =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[faceI],           // modified face
                        faceI,                          // label of face
                        faceOwner[faceI],               // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchI,                         // patch for face
                        false,                          // remove from zone
                        faceZoneI,                      // zone for face
                        flip                            // face flip in zone
                    )
                );
            }
            faceI++;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToZone, cellI)
    {
        label zoneI = cellToZone[cellI];

        if (zoneI >= 0)
        {
            meshMod.setAction
            (
                polyModifyCell
                (
                    cellI,
                    false,          // removeFromZone
                    zoneI
                )
            );
        }
    }
}


void Foam::meshRefinement::allocateInterRegionFaceZone
(
    const label ownZone,
    const label neiZone,
    wordPairHashTable& zonesToFaceZone,
    LabelPairMap<word>& zoneIDsToFaceZone
) const
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    if (ownZone != neiZone)
    {
        // Make sure lowest number cellZone is master. Non-cellZone
        // areas are slave
        const bool swap =
        (
            ownZone == -1
         || (neiZone != -1 && ownZone > neiZone)
        );

        // Quick check whether we already have pair of zones
        labelPair key(ownZone, neiZone);
        if (swap)
        {
            key.flip();
        }

        if (!zoneIDsToFaceZone.found(key))
        {
            // Not found. Allocate.
            const word ownZoneName =
            (
                ownZone != -1
              ? cellZones[ownZone].name()
              : "none"
            );
            const word neiZoneName =
            (
                neiZone != -1
              ? cellZones[neiZone].name()
              : "none"
            );

            // Get lowest zone first
            Pair<word> wordKey(ownZoneName, neiZoneName);
            if (swap)
            {
                wordKey.flip();
            }

            word fzName = wordKey.first() + "_to_" + wordKey.second();

            zoneIDsToFaceZone.insert(key, fzName);
            zonesToFaceZone.insert(wordKey, fzName);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshRefinement::baffleAndSplitMesh
(
    const bool doHandleSnapProblems,
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const label nErodeCellZones,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,

    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh,
    const bool exitIfLeakPath,
    const refPtr<coordSetWriter>& leakPathFormatter,
    refPtr<surfaceWriter>& surfFormatter
)
{
    // Introduce baffles
    // ~~~~~~~~~~~~~~~~~

    // Split the mesh along internal faces wherever there is a pierce between
    // two cell centres.

    Info<< "Introducing baffles for "
        << returnReduce(countHits(), sumOp<label>())
        << " faces that are intersected by the surface." << nl << endl;

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList ownPatch, neiPatch;

    refPtr<surfaceWriter> dummyWriter(nullptr);
    getBafflePatches
    (
        nErodeCellZones,
        globalToMasterPatch,

        locationsInMesh,
        zonesInMesh,
        locationsOutsideMesh,
        exitIfLeakPath,
        refPtr<coordSetWriter>(nullptr),
        dummyWriter,

        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    createBaffles(ownPatch, neiPatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }

    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles", true);

    if (debug&MESH)
    {
        const_cast<Time&>(mesh_.time())++;
        Pout<< "Writing baffled mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"baffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }


    // Introduce baffles to delete problem cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Create some additional baffles where we want surface cells removed.

    if (doHandleSnapProblems)
    {
        handleSnapProblems
        (
            snapParams,
            useTopologicalSnapDetection,
            removeEdgeConnectedCells,
            perpendicularAngle,
            motionDict,
            runTime,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Removing additional cells might have created disconnected bits
        // so re-do the surface intersections
        {
            // Swap neighbouring cell centres and cell level
            neiLevel.setSize(mesh_.nBoundaryFaces());
            neiCc.setSize(mesh_.nBoundaryFaces());
            calcNeighbourData(neiLevel, neiCc);

            labelList ownPatch, neiPatch;
            refPtr<surfaceWriter> dummyWriter(nullptr);
            getBafflePatches
            (
                nErodeCellZones,
                globalToMasterPatch,

                locationsInMesh,
                zonesInMesh,
                locationsOutsideMesh,
                exitIfLeakPath,
                refPtr<coordSetWriter>(nullptr),
                dummyWriter,

                neiLevel,
                neiCc,

                ownPatch,
                neiPatch
            );

            createBaffles(ownPatch, neiPatch);
        }

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            checkData();
        }
    }


    // Select part of mesh
    // ~~~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Remove unreachable sections of mesh" << nl
        << "-----------------------------------" << nl
        << endl;

    if (debug)
    {
        ++runTime;
    }

    splitMeshRegions
    (
        globalToMasterPatch,
        globalToSlavePatch,
        locationsInMesh,
        locationsOutsideMesh,
        true,       // Exit if any connection between inside and outside
        leakPathFormatter
    );

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Split mesh in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After subsetting", true);

    if (debug&MESH)
    {
        ++runTime;
        Pout<< "Writing subsetted mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


void Foam::meshRefinement::mergeFreeStandingBaffles
(
    const bool samePatch,
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const scalar planarAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh
)
{
    // Merge baffles
    // ~~~~~~~~~~~~~

    Info<< nl
        << "Merge free-standing baffles" << nl
        << "---------------------------" << nl
        << endl;


    // List of pairs of freestanding baffle faces.
    List<labelPair> couples
    (
        freeStandingBaffles    // filter out freestanding baffles
        (
            localPointRegion::findDuplicateFacePairs(mesh_),
            planarAngle,
            samePatch
        )
    );

    const label nCouples = returnReduce(couples.size(), sumOp<label>());

    Info<< "Detected free-standing baffles : " << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        mergeBaffles(couples, Map<label>(0));

        // Detect any problem cells resulting from merging of baffles
        // and delete them
        handleSnapProblems
        (
            snapParams,
            useTopologicalSnapDetection,
            removeEdgeConnectedCells,
            perpendicularAngle,
            motionDict,
            runTime,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Very occasionally removing a problem cell might create a disconnected
        // region so re-check

        Info<< nl
            << "Remove unreachable sections of mesh" << nl
            << "-----------------------------------" << nl
            << endl;

        if (debug)
        {
            ++runTime;
        }

        splitMeshRegions
        (
            globalToMasterPatch,
            globalToSlavePatch,
            locationsInMesh,
            locationsOutsideMesh,
            true,   // Exit if any connection between inside and outside
            refPtr<coordSetWriter>(nullptr) // leakPathFormatter
        );


        if (debug)
        {
            // Debug:test all is still synced across proc patches
            checkData();
        }
    }
    Info<< "Merged free-standing baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMesh
(
    const label nBufferLayers,
    const label nErodeCellZones,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,

    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh,
    const bool exitIfLeakPath,
    const refPtr<coordSetWriter>& leakPathFormatter,
    refPtr<surfaceWriter>& surfFormatter
)
{
    // Determine patches to put intersections into
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);

    // Find intersections with all unnamed(!) surfaces
    labelList ownPatch, neiPatch;
    getBafflePatches
    (
        nErodeCellZones,
        globalToMasterPatch,

        locationsInMesh,
        zonesInMesh,
        locationsOutsideMesh,
        exitIfLeakPath,
        leakPathFormatter,
        surfFormatter,

        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces(), false);

    forAll(ownPatch, faceI)
    {
        if (ownPatch[faceI] != -1 || neiPatch[faceI] != -1)
        {
            blockedFace[faceI] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());


    regionSplit cellRegion(mesh_, blockedFace);

    // Set unreachable cells to -1
    findRegions
    (
        mesh_,
        vector::uniform(mergeDistance_),   // perturbVec
        locationsInMesh,
        locationsOutsideMesh,
        cellRegion.nRegions(),
        cellRegion,
        blockedFace,
        // Leak-path
        false,      // do not exit if outside location found
        leakPathFormatter
    );

    return splitMesh
    (
        nBufferLayers,
        globalToMasterPatch,
        globalToSlavePatch,
        cellRegion,
        ownPatch,
        neiPatch
    );
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMesh
(
    const label nBufferLayers,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,

    labelList& cellRegion,
    labelList& ownPatch,
    labelList& neiPatch
)
{
    // Walk out nBufferlayers from region boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (modifies cellRegion, ownPatch)
    // Takes over face patch onto points and then back to faces and cells
    // (so cell-face-point walk)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    // Checks
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (ownPatch[facei] == -1 && neiPatch[facei] != -1)
        {
            FatalErrorInFunction << "Problem in face:" << facei
                << " at:" << mesh_.faceCentres()[facei]
                << " ownPatch:" << ownPatch[facei]
                << " neiPatch:" << neiPatch[facei]
                << exit(FatalError);
        }
        else
        {
            // Check if cellRegion indeed limited by patch
            const label ownRegion = cellRegion[faceOwner[facei]];
            const label neiRegion = cellRegion[faceNeighbour[facei]];
            if (ownRegion != neiRegion)
            {
                if (ownPatch[facei] == -1)
                {
                    FatalErrorInFunction << "Problem in face:" << facei
                        << " at:" << mesh_.faceCentres()[facei]
                        << " ownPatch:" << ownPatch[facei]
                        << " ownRegion:" << ownRegion
                        << " neiPatch:" << neiPatch[facei]
                        << " neiRegion:" << neiRegion
                        << exit(FatalError);
                }
            }
        }
    }

    // Determine on original data the nearest face. This is used as a fall-back
    // to set the patch if the nBufferLayers walking didn't work.
    labelList nearestOwnPatch;
    if (nBufferLayers)
    {
        DynamicList<label> startFaces;
        forAll(ownPatch, facei)
        {
            if (ownPatch[facei] != -1)
            {
                startFaces.append(facei);
            }
        }

        // Per face the index to the start face.
        labelList faceToStart;
        autoPtr<mapDistribute> mapPtr;
        nearestFace
        (
            startFaces,
            bitSet(mesh_.nFaces()), // no blocked faces
            mapPtr,
            faceToStart,
            nBufferLayers+4     // bit more than nBufferLayers since
                                // walking face-cell-face
        );

        // Use map to push ownPatch to all reached faces
        labelList startOwnPatch(ownPatch, startFaces);
        mapPtr().distribute(startOwnPatch);

        nearestOwnPatch.setSize(mesh_.nFaces());
        nearestOwnPatch = -1;
        forAll(faceToStart, facei)
        {
            const label sloti = faceToStart[facei];
            if (sloti != -1)
            {
                nearestOwnPatch[facei] = startOwnPatch[sloti];
            }
        }
    }

    // Leak closure:
    // ~~~~~~~~~~~~~
    // We do not want to add buffer layers on the frozen points/faces
    // since these are the exact faces needed to close a hole (to an
    // locationOutsideMesh). Adding even
    // a single layer of cells would mean that in further manipulation there
    // is now no path to the locationOutsideMesh so the layer closure does
    // not get triggered and we keep the added 1 layer of cells on the
    // closure faces.
    bitSet isFrozenPoint(mesh_.nPoints());
    bitSet isFrozenFace(mesh_.nFaces());

    if (nBufferLayers)
    {
        const labelListList& pointFaces = mesh_.pointFaces();
        const pointZoneMesh& pzs = mesh_.pointZones();
        const label pointZonei = pzs.findZoneID("frozenPoints");
        if (pointZonei != -1)
        {
            const pointZone& pz = pzs[pointZonei];
            isFrozenPoint.set(pz);
            for (const label pointi : pz)
            {
                isFrozenFace.set(pointFaces[pointi]);
            }
        }

        const faceZoneMesh& fzs = mesh_.faceZones();
        const label faceZonei = fzs.findZoneID("frozenFaces");
        if (faceZonei != -1)
        {
            const faceZone& fz = fzs[faceZonei];
            isFrozenFace.set(fz);
            for (const label facei : fz)
            {
                isFrozenPoint.set(mesh_.faces()[facei]);
            }
        }
    }

    // Patch for exposed faces for lack of anything sensible.
    label defaultPatch = 0;
    if (globalToMasterPatch.size())
    {
        defaultPatch = globalToMasterPatch[0];
    }

    for (label i = 0; i < nBufferLayers; i++)
    {
        // 1. From cells (via faces) to points

        labelList pointBaffle(mesh_.nPoints(), -1);

        forAll(faceNeighbour, facei)
        {
            if (!isFrozenFace[facei])
            {
                const face& f = mesh_.faces()[facei];

                const label ownRegion = cellRegion[faceOwner[facei]];
                const label neiRegion = cellRegion[faceNeighbour[facei]];

                if (ownRegion == -1 && neiRegion != -1)
                {
                    // Note max(..) since possibly regionSplit might have split
                    // off extra unreachable parts of mesh. Note: or can this
                    // only happen for boundary faces?
                    forAll(f, fp)
                    {
                        if (!isFrozenPoint[f[fp]])
                        {
                            pointBaffle[f[fp]] =
                                max(defaultPatch, ownPatch[facei]);
                        }
                    }
                }
                else if (ownRegion != -1 && neiRegion == -1)
                {
                    label newPatchi = neiPatch[facei];
                    if (newPatchi == -1)
                    {
                        newPatchi = max(defaultPatch, ownPatch[facei]);
                    }
                    forAll(f, fp)
                    {
                        if (!isFrozenPoint[f[fp]])
                        {
                            pointBaffle[f[fp]] = newPatchi;
                        }
                    }
                }
            }
        }

        labelList neiCellRegion;
        syncTools::swapBoundaryCellList(mesh_, cellRegion, neiCellRegion);

        //for
        //(
        //    label facei = mesh_.nInternalFaces();
        //    facei < mesh_.nFaces();
        //    facei++
        //)
        //{
        //    if (!isFrozenFace[facei])
        //    {
        //        const face& f = mesh_.faces()[facei];
        //
        //        const label ownRegion = cellRegion[faceOwner[facei]];
        //        const label neiRegion =
        //            neiCellRegion[facei-mesh_.nInternalFaces()];
        //
        //        if (ownRegion == -1 && neiRegion != -1)
        //        {
        //            forAll(f, fp)
        //            {
        //                if (!isFrozenPoint[f[fp]])
        //                {
        //                    pointBaffle[f[fp]] =
        //                        max(defaultPatch, ownPatch[facei]);
        //                }
        //            }
        //        }
        //    }
        //}

        const auto& patches = mesh_.boundaryMesh();
        for (const auto& pp : patches)
        {
            if (pp.coupled())
            {
                // Note: swapBoundaryCellList only works on cyclic&processor.
                //       Does not handle e.g. cyclicAMI. TBD?
                // Note: we could check check our side being in the set
                //       since syncPointList below will push over any decision
                //       made by the other side.
                forAll(pp, i)
                {
                    const label facei = pp.start()+i;
                    if (!isFrozenFace[facei])
                    {
                        const face& f = mesh_.faces()[facei];
                        const label ownRegion = cellRegion[faceOwner[facei]];
                        const label neiRegion =
                            neiCellRegion[facei-mesh_.nInternalFaces()];

                        // Same logic as for internal faces

                        if (ownRegion == -1 && neiRegion != -1)
                        {
                            forAll(f, fp)
                            {
                                if (!isFrozenPoint[f[fp]])
                                {
                                    pointBaffle[f[fp]] =
                                        max(defaultPatch, ownPatch[facei]);
                                }
                            }
                        }
                        else if (ownRegion != -1 && neiRegion == -1)
                        {
                            label newPatchI = neiPatch[facei];
                            if (newPatchI == -1)
                            {
                                newPatchI = max(defaultPatch, ownPatch[facei]);
                            }
                            forAll(f, fp)
                            {
                                if (!isFrozenPoint[f[fp]])
                                {
                                    pointBaffle[f[fp]] = newPatchI;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                forAll(pp, i)
                {
                    const label facei = pp.start()+i;
                    if (!isFrozenFace[facei])
                    {
                        const face& f = mesh_.faces()[facei];
                        const label ownRegion = cellRegion[faceOwner[facei]];

                        if (ownRegion != -1)
                        {
                            forAll(f, fp)
                            {
                                if (!isFrozenPoint[f[fp]])
                                {
                                    pointBaffle[f[fp]] =
                                        max(defaultPatch, ownPatch[facei]);
                                }
                            }
                        }
                    }
                }
            }
        }


        // Sync
        syncTools::syncPointList
        (
            mesh_,
            pointBaffle,
            maxEqOp<label>(),
            label(-1)           // null value
        );


        // 2. From points back to faces

        const labelListList& pointFaces = mesh_.pointFaces();

        forAll(pointFaces, pointi)
        {
            if (pointBaffle[pointi] != -1)
            {
                const labelList& pFaces = pointFaces[pointi];

                forAll(pFaces, pFacei)
                {
                    const label facei = pFaces[pFacei];

                    if (!isFrozenFace[facei] && ownPatch[facei] == -1)
                    {
                        ownPatch[facei] = pointBaffle[pointi];
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());


        // 3. From faces to cells (cellRegion) and back to faces (ownPatch)

        labelList newOwnPatch(ownPatch);

        forAll(ownPatch, facei)
        {
            if (!isFrozenFace[facei] && ownPatch[facei] != -1)
            {
                const label own = faceOwner[facei];

                if (cellRegion[own] == -1)
                {
                    cellRegion[own] = labelMax;

                    const cell& ownFaces = mesh_.cells()[own];
                    forAll(ownFaces, j)
                    {
                        const label ownFacei = ownFaces[j];
                        if (!isFrozenFace[ownFacei] && ownPatch[ownFacei] == -1)
                        {
                            newOwnPatch[ownFacei] = ownPatch[facei];
                        }
                    }
                }
                if (mesh_.isInternalFace(facei))
                {
                    const label nei = faceNeighbour[facei];

                    if (cellRegion[nei] == -1)
                    {
                        cellRegion[nei] = labelMax;

                        const cell& neiFaces = mesh_.cells()[nei];
                        forAll(neiFaces, j)
                        {
                            const label neiFacei = neiFaces[j];
                            const bool isFrozen = isFrozenFace[neiFacei];
                            if (!isFrozen && ownPatch[neiFacei] == -1)
                            {
                                newOwnPatch[neiFacei] = ownPatch[facei];
                            }
                        }
                    }
                }
            }
        }

        ownPatch.transfer(newOwnPatch);

        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    }


    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] == -1)
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    const label nCellsToKeep = returnReduce
    (
        mesh_.nCells() - cellsToRemove.size(),
        sumOp<label>()
    );

    Info<< "Keeping all cells containing inside points" << endl
        << "Selected for keeping : " << nCellsToKeep << " cells." << endl;


    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    const labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatches(exposedFaces.size());

    label nUnpatched = 0;

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        if (ownPatch[facei] != -1)
        {
            exposedPatches[i] = ownPatch[facei];
        }
        else
        {
            const label fallbackPatch =
            (
                nearestOwnPatch.size()
              ? nearestOwnPatch[facei]
              : defaultPatch
            );
            if (nUnpatched == 0)
            {
                WarningInFunction
                    << "For exposed face " << facei
                    << " fc:" << mesh_.faceCentres()[facei]
                    << " found no patch." << endl
                    << "    Taking patch " << fallbackPatch
                    << " instead. Suppressing future warnings" << endl;
            }
            nUnpatched++;

            exposedPatches[i] = fallbackPatch;
        }
    }

    reduce(nUnpatched, sumOp<label>());
    if (nUnpatched > 0)
    {
        Info<< "Detected " << nUnpatched << " faces out of "
            << returnReduce(exposedFaces.size(), sumOp<label>())
            << " for which the default patch " << defaultPatch
            << " will be used" << endl;
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::removeLimitShells
(
    const label nBufferLayers,
    const label nErodeCellZones,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh
)
{
    // Determine patches to put intersections into
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);

    // Find intersections with all unnamed(!) surfaces
    labelList ownPatch, neiPatch;
    refPtr<surfaceWriter> dummyWriter(nullptr);
    getBafflePatches
    (
        nErodeCellZones,
        globalToMasterPatch,

        locationsInMesh,
        zonesInMesh,
        locationsOutsideMesh,
        false,                      // do not exit. Use leak-closure instead.
        refPtr<coordSetWriter>(nullptr),
        dummyWriter,

        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );


    labelList cellRegion(mesh_.nCells(), Zero);
    // Find any cells inside a limit shell with minLevel -1
    labelList levelShell;
    limitShells_.findLevel
    (
        mesh_.cellCentres(),
        labelList(mesh_.nCells(), -1),  // pick up only shells with -1
        levelShell
    );
    forAll(levelShell, celli)
    {
        if (levelShell[celli] != -1)
        {
            // Mark cell region so it gets deleted
            cellRegion[celli] = -1;
        }
    }

    autoPtr<mapPolyMesh> mapPtr = splitMesh
    (
        nBufferLayers,
        globalToMasterPatch,
        globalToSlavePatch,
        cellRegion,
        ownPatch,
        neiPatch
    );

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh_.time())++;
        Pout<< "Writing mesh after removing limitShells"
            << " to time " << timeName() << endl;
        write
        (
            debugType(debug),
            writeType
            (
                writeLevel()
              | WRITEMESH
            ),
            mesh_.time().path()/timeName()
        );
    }
    return mapPtr;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints
(
    const localPointRegion& regionSide
)
{
    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nNonManifPoints = returnReduce
    (
        regionSide.meshPointMap().size(),
        sumOp<label>()
    );

    Info<< "dupNonManifoldPoints : Found : " << nNonManifPoints
        << " non-manifold points (out of "
        << mesh_.globalData().nTotalPoints()
        << ')' << endl;


    autoPtr<mapPolyMesh> mapPtr;

    if (nNonManifPoints)
    {
        // Topo change engine
        duplicatePoints pointDuplicator(mesh_);

        // Insert changes into meshMod
        pointDuplicator.setRefinement(regionSide, meshMod);

        // Remove any unnecessary fields
        mesh_.clearOut();
        mesh_.moving(false);

        // Change the mesh (no inflation, parallel sync)
        mapPtr = meshMod.changeMesh(mesh_, false, true);
        mapPolyMesh& map = *mapPtr;

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map.hasMotionPoints())
        {
            mesh_.movePoints(map.preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Is mapping only (no faces created, positions
        // stay same) so no need to recalculate intersections.
        updateMesh(map, labelList());
    }

    return mapPtr;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints()
{
    // Analyse which points need to be duplicated
    localPointRegion regionSide(mesh_);

    return dupNonManifoldPoints(regionSide);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergePoints
(
    const labelList& pointToDuplicate
)
{
    label nPointPairs = 0;
    forAll(pointToDuplicate, pointI)
    {
        label otherPointI = pointToDuplicate[pointI];
        if (otherPointI != -1)
        {
            nPointPairs++;
        }
    }

    autoPtr<mapPolyMesh> mapPtr;

    if (returnReduceOr(nPointPairs))
    {
        Map<label> pointToMaster(2*nPointPairs);
        forAll(pointToDuplicate, pointI)
        {
            label otherPointI = pointToDuplicate[pointI];
            if (otherPointI != -1)
            {
                // Slave point
                pointToMaster.insert(pointI, otherPointI);
            }
        }

        // Topochange container
        polyTopoChange meshMod(mesh_);

        // Insert changes
        polyMeshAdder::mergePoints(mesh_, pointToMaster, meshMod);

        // Remove any unnecessary fields
        mesh_.clearOut();
        mesh_.moving(false);

        // Change the mesh (no inflation, parallel sync)
        mapPtr = meshMod.changeMesh(mesh_, false, true);
        mapPolyMesh& map = *mapPtr;

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map.hasMotionPoints())
        {
            mesh_.movePoints(map.preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Is mapping only (no faces created, positions
        // stay same) so no need to recalculate intersections.
        updateMesh(map, labelList());
    }

    return mapPtr;
}


// Duplicate points on 'boundary' zones. Do not duplicate points on
// 'internal' or 'baffle' zone. Whether points are on normal patches does
// not matter
Foam::autoPtr<Foam::mapPolyMesh>
Foam::meshRefinement::dupNonManifoldBoundaryPoints()
{
    const labelList boundaryFaceZones
    (
        getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::BOUNDARY
            )
        )
    );
    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = getZones(fzTypes);
    }



    // 0 : point used by normal, unzoned boundary faces
    // 1 : point used by 'boundary' zone
    // 2 : point used by internal/baffle zone
    PackedList<2> pointStatus(mesh_.nPoints(), 0u);

    forAll(boundaryFaceZones, j)
    {
        const faceZone& fZone = mesh_.faceZones()[boundaryFaceZones[j]];
        forAll(fZone, i)
        {
            const face& f = mesh_.faces()[fZone[i]];
            forAll(f, fp)
            {
                pointStatus[f[fp]] = max(pointStatus[f[fp]], 1u);
            }
        }
    }
    forAll(internalOrBaffleFaceZones, j)
    {
        const faceZone& fZone = mesh_.faceZones()[internalOrBaffleFaceZones[j]];
        forAll(fZone, i)
        {
            const face& f = mesh_.faces()[fZone[i]];
            forAll(f, fp)
            {
                pointStatus[f[fp]] = max(pointStatus[f[fp]], 2u);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pointStatus,
        maxEqOp<unsigned int>(),    // combine op
        0u                          // null value
    );

    // Pick up points on boundary zones that are not on internal/baffle zones
    label n = 0;
    forAll(pointStatus, pointI)
    {
        if (pointStatus[pointI] == 1u)
        {
            n++;
        }
    }

    label globalNPoints = returnReduce(n, sumOp<label>());
    Info<< "Duplicating " << globalNPoints << " points on"
        << " faceZones of type "
        << surfaceZonesInfo::faceZoneTypeNames[surfaceZonesInfo::BOUNDARY]
        << endl;

    autoPtr<mapPolyMesh> map;

    if (globalNPoints)
    {
        labelList candidatePoints(n);
        n = 0;
        forAll(pointStatus, pointI)
        {
            if (pointStatus[pointI] == 1u)
            {
                candidatePoints[n++] = pointI;
            }
        }
        localPointRegion regionSide(mesh_, candidatePoints);
        map = dupNonManifoldPoints(regionSide);
    }
    return map;
}


// Zoning
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::zonify
(
    const bool allowFreeStandingZoneFaces,
    const label nErodeCellZones,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh,
    const bool exitIfLeakPath,
    const refPtr<coordSetWriter>& leakPathFormatter,
    refPtr<surfaceWriter>& surfFormatter,
    wordPairHashTable& zonesToFaceZone
)
{
    if (locationsInMesh.size() != zonesInMesh.size())
    {
        FatalErrorInFunction << "problem" << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();


    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);


    // Add any faceZones, cellZones originating from surface to the mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList surfaceToCellZone;
    labelListList surfaceToFaceZones;

    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));
    if (namedSurfaces.size())
    {
        Info<< "Setting cellZones according to named surfaces:" << endl;
        forAll(namedSurfaces, i)
        {
            label surfI = namedSurfaces[i];

            Info<< "Surface : " << surfaces_.names()[surfI] << nl
                << "    faceZones : " << surfZones[surfI].faceZoneNames() << nl
                << "    cellZone  : " << surfZones[surfI].cellZoneName()
                << endl;
        }
        Info<< endl;

        // Add zones to mesh
        surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
        surfaceToFaceZones = surfaceZonesInfo::addFaceZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
    }



    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Zone per cell:
    // -2 : unset : not allowed!
    // -1 : not in any zone (zone 'none')
    // >=0: zoneID
    // namedSurfaceRegion: zero sized or:
    // -1  : face not intersecting a named surface
    // >=0 : index of named surface
    labelList cellToZone;
    labelList namedSurfaceRegion;
    bitSet posOrientation;
    {
        labelList unnamedRegion1;
        labelList unnamedRegion2;

        zonify
        (
            allowFreeStandingZoneFaces,
            nErodeCellZones,// Use erosion (>0) or regionSplit to clean up
            -1,             // Set all cells with cellToZone -2 to -1
            locationsInMesh,
            zonesInMesh,
            locationsOutsideMesh,
            exitIfLeakPath,
            leakPathFormatter,
            surfFormatter,

            cellToZone,
            unnamedRegion1,
            unnamedRegion2,
            namedSurfaceRegion,
            posOrientation
        );
    }


    // Convert namedSurfaceRegion (index of named surfaces) to
    // actual faceZone index

    //- Per face index of faceZone or -1
    labelList faceToZone(mesh_.nFaces(), -1);

    forAll(namedSurfaceRegion, faceI)
    {
        //label surfI = namedSurfaceIndex[faceI];
        label globalI = namedSurfaceRegion[faceI];
        if (globalI != -1)
        {
            const labelPair spr = surfaces_.whichSurface(globalI);
            faceToZone[faceI] = surfaceToFaceZones[spr.first()][spr.second()];
        }
    }



    // Allocate and assign faceZones from cellZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // 1. Detect inter-region face and allocate names

        LabelPairMap<word> zoneIDsToFaceZone;

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            if (faceToZone[faceI] == -1)    // Not named surface
            {
                // Face not yet in a faceZone. (it might already have been
                // done so by a 'named' surface). Check if inbetween different
                // cellZones
                allocateInterRegionFaceZone
                (
                    cellToZone[mesh_.faceOwner()[faceI]],
                    cellToZone[mesh_.faceNeighbour()[faceI]],
                    zonesToFaceZone,
                    zoneIDsToFaceZone
                );
            }
        }

        labelList neiCellZone;
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

        forAll(neiCellZone, bFaceI)
        {
            label faceI = bFaceI + mesh_.nInternalFaces();
            if (faceToZone[faceI] == -1)
            {
                allocateInterRegionFaceZone
                (
                    cellToZone[mesh_.faceOwner()[faceI]],
                    neiCellZone[bFaceI],
                    zonesToFaceZone,
                    zoneIDsToFaceZone
                );
            }
        }


        // 2.Combine faceZoneNames allocated on different processors

        Pstream::mapCombineReduce(zonesToFaceZone, eqOp<word>());


        // 3. Allocate faceZones from (now synchronised) faceZoneNames
        //    Note: the faceZoneNames contain the same data but in different
        //          order. We could sort the contents but instead just loop
        //          in sortedToc order.

        Info<< "Setting faceZones according to neighbouring cellZones:"
            << endl;

        // From cellZone indices to faceZone index
        LabelPairMap<label> fZoneLookup(zonesToFaceZone.size());

        const cellZoneMesh& cellZones = mesh_.cellZones();

        {
            List<Pair<word>> czs(zonesToFaceZone.sortedToc());

            forAll(czs, i)
            {
                const Pair<word>& cz = czs[i];
                const word& fzName = zonesToFaceZone[cz];

                Info<< indent<< "cellZones : "
                    << cz[0] << ' ' << cz[1] << nl
                    << "    faceZone : " << fzName << endl;

                label faceZoneI = surfaceZonesInfo::addFaceZone
                (
                    fzName,                 // name
                    labelList(0),           // addressing
                    boolList(0),            // flipMap
                    mesh_
                );

                label cz0 = cellZones.findZoneID(cz[0]);
                label cz1 = cellZones.findZoneID(cz[1]);

                fZoneLookup.insert(labelPair(cz0, cz1), faceZoneI);
            }
        }


        // 4. Set faceToZone with new faceZones


        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            if (faceToZone[faceI] == -1)
            {
                // Face not yet in a faceZone. (it might already have been
                // done so by a 'named' surface). Check if inbetween different
                // cellZones

                label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
                label neiZone = cellToZone[mesh_.faceNeighbour()[faceI]];
                if (ownZone != neiZone)
                {
                    const bool swap =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                    labelPair key(ownZone, neiZone);
                    if (swap)
                    {
                        key.flip();
                    }
                    faceToZone[faceI] = fZoneLookup[key];
                }
            }
        }
        forAll(neiCellZone, bFaceI)
        {
            label faceI = bFaceI + mesh_.nInternalFaces();
            if (faceToZone[faceI] == -1)
            {
                label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
                label neiZone = neiCellZone[bFaceI];
                if (ownZone != neiZone)
                {
                    const bool swap =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                    labelPair key(ownZone, neiZone);
                    if (swap)
                    {
                        key.flip();
                    }
                    faceToZone[faceI] = fZoneLookup[key];
                }
            }
        }
        Info<< endl;
    }




    // Get coupled neighbour cellZone. Set to -1 on non-coupled patches.
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!pp.coupled())
        {
            label bFaceI = pp.start()-mesh_.nInternalFaces();
            forAll(pp, i)
            {
                neiCellZone[bFaceI++] = -1;
            }
        }
    }



    // Get per face whether is it master (of a coupled set of faces)
    const bitSet isMasterFace(syncTools::getMasterFaces(mesh_));


    // faceZones
    // ~~~~~~~~~
    // Faces on faceZones come in two variants:
    // - faces on the outside of a cellZone. They will be oriented to
    //   point out of the maximum cellZone.
    // - free-standing faces. These will be oriented according to the
    //   local surface normal. We do this in a two step algorithm:
    //      - do a consistent orientation
    //      - check number of faces with consistent orientation
    //      - if <0 flip the whole patch
    bitSet meshFlipMap(mesh_.nFaces(), false);
    {
        // Collect all data on zone faces without cellZones on either side.
        const indirectPrimitivePatch patch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                freeStandingBaffleFaces
                (
                    faceToZone,
                    cellToZone,
                    neiCellZone
                )
            ),
            mesh_.points()
        );

        label nFreeStanding = returnReduce(patch.size(), sumOp<label>());
        if (nFreeStanding > 0)
        {
            Info<< "Detected " << nFreeStanding << " free-standing zone faces"
                << endl;

            if (debug)
            {
                OBJstream str(mesh_.time().path()/"freeStanding.obj");
                Pout<< "meshRefinement::zonify : dumping faceZone faces to "
                    << str.name() << endl;
                str.write(patch.localFaces(), patch.localPoints(), false);
            }


            // Detect non-manifold edges
            labelList nMasterFacesPerEdge;
            calcPatchNumMasterFaces(isMasterFace, patch, nMasterFacesPerEdge);

            // Mark zones. Even a single original surface might create multiple
            // disconnected/non-manifold-connected zones
            labelList faceToConnectedZone;
            const label nZones = markPatchZones
            (
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone
            );

            Map<label> nPosOrientation(2*nZones);
            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                nPosOrientation.insert(zoneI, 0);
            }

            // Make orientations consistent in a topological way. This just
            // checks  the first face per zone for whether nPosOrientation
            // is negative (which it never is at this point)
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );

            // Count per region the number of orientations (taking the new
            // flipMap into account)
            forAll(patch.addressing(), faceI)
            {
                label meshFaceI = patch.addressing()[faceI];

                if (isMasterFace.test(meshFaceI))
                {
                    label n = 1;
                    if
                    (
                        posOrientation.test(meshFaceI)
                     == meshFlipMap.test(meshFaceI)
                    )
                    {
                        n = -1;
                    }

                    nPosOrientation.find(faceToConnectedZone[faceI])() += n;
                }
            }
            Pstream::mapCombineReduce(nPosOrientation, plusEqOp<label>());


            Info<< "Split " << nFreeStanding << " free-standing zone faces"
                << " into " << nZones << " disconnected regions with size"
                << " (negative denotes wrong orientation) :"
                << endl;

            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                Info<< "    " << zoneI << "\t" << nPosOrientation[zoneI]
                    << endl;
            }
            Info<< endl;


            // Re-apply with new counts (in nPosOrientation). This will cause
            // zones with a negative count to be flipped.
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );
        }
    }




    // Topochange container
    polyTopoChange meshMod(mesh_);

    // Insert changes to put cells and faces into zone
    zonify
    (
        isMasterFace,
        cellToZone,
        neiCellZone,
        faceToZone,
        meshFlipMap,
        meshMod
    );

    // Remove any unnecessary fields
    mesh_.clearOut();
    mesh_.moving(false);

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map());

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Print some stats (note: zones are synchronised)
    if (mesh_.cellZones().size() > 0)
    {
        Info<< "CellZones:" << endl;
        forAll(mesh_.cellZones(), zoneI)
        {
            const cellZone& cz = mesh_.cellZones()[zoneI];
            Info<< "    " << cz.name()
                << "\tsize:" << returnReduce(cz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }
    if (mesh_.faceZones().size() > 0)
    {
        Info<< "FaceZones:" << endl;
        forAll(mesh_.faceZones(), zoneI)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            Info<< "    " << fz.name()
                << "\tsize:" << returnReduce(fz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }

    // None of the faces has changed, only the zones. Still...
    updateMesh(map(), labelList());

    return map;
}


// ************************************************************************* //
