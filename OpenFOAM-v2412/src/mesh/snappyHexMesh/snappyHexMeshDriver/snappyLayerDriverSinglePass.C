/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Description
    Single pass layer addition. Can be removed once multi-pass works ok.

\*----------------------------------------------------------------------------*/

#include "snappyLayerDriver.H"
//#include "motionSmoother.H"
//#include "pointSet.H"
//#include "faceSet.H"
//#include "cellSet.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "addPatchCellLayer.H"
#include "mapDistributePolyMesh.H"
//#include "OBJstream.H"
#include "layerParameters.H"
#include "externalDisplacementMeshMover.H"
//#include "profiling.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::snappyLayerDriver::addLayersSinglePass
(
    const layerParameters& layerParams,
    const dictionary& motionDict,
    const labelList& patchIDs,
    const label nAllowableErrors,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();

    // faceZones of type internal or baffle (for merging points across)
    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = meshRefiner_.getZones(fzTypes);
    }

    // faceZones of type internal (for checking mesh quality across and
    // merging baffles)
    const labelList internalFaceZones
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::INTERNAL
            )
        )
    );

    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    List<labelPair> baffles;
    {
        labelList originatingFaceZone;
        meshRefiner_.createZoneBaffles
        (
            identity(mesh.faceZones().size()),
            baffles,
            originatingFaceZone
        );

        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing baffled mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName()
            );
        }
    }


    // Duplicate points on faceZones of type boundary. Should normally already
    // be done by snapping phase
    {
        autoPtr<mapPolyMesh> map = meshRefiner_.dupNonManifoldBoundaryPoints();
        if (map)
        {
            const labelList& reverseFaceMap = map->reverseFaceMap();
            forAll(baffles, i)
            {
                label f0 = reverseFaceMap[baffles[i].first()];
                label f1 = reverseFaceMap[baffles[i].second()];
                baffles[i] = labelPair(f0, f1);
            }
        }
    }



    // Now we have all patches determine the number of layer per patch
    // This will be layerParams.numLayers() except for faceZones where one
    // side does get layers and the other not in which case we want to
    // suppress movement by explicitly setting numLayers 0
    labelList numLayers(layerParams.numLayers());
    {
        labelHashSet layerIDs(patchIDs);
        forAll(mesh.faceZones(), zonei)
        {
            label mpi, spi;
            surfaceZonesInfo::faceZoneType fzType;
            bool hasInfo = meshRefiner_.getFaceZoneInfo
            (
                mesh.faceZones()[zonei].name(),
                mpi,
                spi,
                fzType
            );
            if (hasInfo)
            {
                const polyBoundaryMesh& pbm = mesh.boundaryMesh();
                if (layerIDs.found(mpi) && !layerIDs.found(spi))
                {
                    // Only layers on master side. Fix points on slave side
                    Info<< "On faceZone " << mesh.faceZones()[zonei].name()
                        << " adding layers to master patch " << pbm[mpi].name()
                        << " only. Freezing points on slave patch "
                        << pbm[spi].name() << endl;
                    numLayers[spi] = 0;
                }
                else if (!layerIDs.found(mpi) && layerIDs.found(spi))
                {
                    // Only layers on slave side. Fix points on master side
                    Info<< "On faceZone " << mesh.faceZones()[zonei].name()
                        << " adding layers to slave patch " << pbm[spi].name()
                        << " only. Freezing points on master patch "
                        << pbm[mpi].name() << endl;
                    numLayers[mpi] = 0;
                }
            }
        }
    }


    // Duplicate points on faceZones that layers are added to
    labelList pointToMaster;
    dupFaceZonePoints
    (
        meshRefiner_,
        patchIDs,  // patch indices
        numLayers, // number of layers per patch
        baffles,
        pointToMaster
    );


    // Add layers to patches
    // ~~~~~~~~~~~~~~~~~~~~~

    // Now we have
    // - mesh with optional baffles and duplicated points for faceZones
    //   where layers are to be added
    // - pointToMaster    : correspondence for duplicated points
    // - baffles          : list of pairs of faces


    autoPtr<indirectPrimitivePatch> pp
    (
        meshRefinement::makePatch
        (
            mesh,
            patchIDs
        )
    );


    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches). This is used only to string up edges between coupled
    // faces (all edges between same (global)face indices get extruded).
    const labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            *pp
        )
    );

    // Determine patches for extruded boundary edges. Calculates if any
    // additional processor patches need to be constructed.

    labelList edgePatchID;
    labelList edgeZoneID;
    boolList edgeFlip;
    labelList inflateFaceID;
    determineSidePatches
    (
        meshRefiner_,
        globalFaces,
        edgeGlobalFaces,
        *pp,

        edgePatchID,
        edgeZoneID,
        edgeFlip,
        inflateFaceID
    );


    // Point-wise extrusion data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Displacement for all pp.localPoints. Set to large value so does
    // not trigger the minThickness truncation (see syncPatchDisplacement below)
    vectorField patchDisp(pp().nPoints(), vector(GREAT, GREAT, GREAT));

    // Number of layers for all pp.localPoints. Note: only valid if
    // extrudeStatus = EXTRUDE.
    labelList patchNLayers(pp().nPoints(), Zero);

    // Whether to add edge for all pp.localPoints.
    List<extrudeMode> extrudeStatus(pp().nPoints(), EXTRUDE);

    // Ideal number of cells added
    const label nIdealTotAddedCells = setPointNumLayers
    (
        layerParams,

        numLayers,
        patchIDs,
        pp(),
        edgeGlobalFaces,

        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    // Determine (wanted) point-wise overall layer thickness and expansion
    // ratio
    scalarField thickness(pp().nPoints());
    scalarIOField minThickness
    (
        IOobject
        (
            "minThickness",
            meshRefiner_.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        pp().nPoints()
    );
    scalarField expansionRatio(pp().nPoints());
    calculateLayerThickness
    (
        *pp,
        patchIDs,
        layerParams,
        meshRefiner_.meshCutter().cellLevel(),
        patchNLayers,
        edge0Len,

        thickness,
        minThickness,
        expansionRatio
    );



    // Per cell 0 or number of layers in the cell column it is part of
    labelList cellNLayers;
    // Per face actual overall layer thickness
    scalarField faceRealThickness;
    // Per face wanted overall layer thickness
    scalarField faceWantedThickness(mesh.nFaces(), Zero);
    {
        UIndirectList<scalar>(faceWantedThickness, pp->addressing()) =
            avgPointData(*pp, thickness);
    }


    // Current set of topology changes. (changing mesh clears out
    // polyTopoChange)
    polyTopoChange meshMod(mesh.boundaryMesh().size());

    // Play changes into meshMod
    addLayers
    (
        layerParams,
        layerParams.nLayerIter(),

        // Mesh quality
        motionDict,
        layerParams.nRelaxedIter(),
        nAllowableErrors,

        patchIDs,
        internalFaceZones,
        baffles,
        numLayers,
        nIdealTotAddedCells,

        // Patch information
        globalFaces,
        pp(),
        edgeGlobalFaces,
        edgePatchID,
        edgeZoneID,
        edgeFlip,
        inflateFaceID,

        // Per patch point the wanted thickness
        thickness,
        minThickness,
        expansionRatio,

        // Per patch point the obtained thickness
        patchDisp,
        patchNLayers,
        extrudeStatus,

        // Complete mesh changes
        meshMod,

        // Stats
        cellNLayers,
        faceRealThickness
    );


    // At this point we have a (shrunk) mesh and a set of topology changes
    // which will make a valid mesh with layer. Apply these changes to the
    // current mesh.

    {
        // Apply the stored topo changes to the current mesh.
        autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);

        mapPolyMesh& map = *mapPtr;

        // Hack to remove meshPhi - mapped incorrectly. TBD.
        mesh.moving(false);
        mesh.clearOut();

        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map.hasMotionPoints())
        {
            mesh.movePoints(map.preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh.setInstance(meshRefiner_.timeName());

        meshRefiner_.updateMesh(map, labelList(0));

        // Update numbering of faceWantedThickness
        meshRefinement::updateList
        (
            map.faceMap(),
            scalar(0),
            faceWantedThickness
        );

        // Print data now that we still have patches for the zones
        //if (meshRefinement::outputLevel() & meshRefinement::OUTPUTLAYERINFO)
        printLayerData
        (
            mesh,
            patchIDs,
            cellNLayers,
            faceWantedThickness,
            faceRealThickness,
            layerParams
        );


        // Dump for debugging
        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing mesh with layers but disconnected to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName()
            );
        }


        // Map baffles, pointToMaster
        mapFaceZonePoints(meshRefiner_, map, baffles, pointToMaster);
    }


    // Merge baffles
    mergeFaceZonePoints
    (
        pointToMaster, // -1 or index of duplicated point
        cellNLayers,
        faceRealThickness,
        faceWantedThickness
    );


    // Do final balancing
    // ~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Balance. No restriction on face zones and baffles.
        autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
        (
            false,
            false,
            labelList::null(),
            scalarField(mesh.nCells(), 1.0),
            decomposer,
            distributor
        );

        // Re-distribute flag of layer faces and cells
        map().distributeCellData(cellNLayers);
        map().distributeFaceData(faceWantedThickness);
        map().distributeFaceData(faceRealThickness);
    }


    // Write mesh data
    // ~~~~~~~~~~~~~~~

    if (!dryRun_)
    {
        writeLayerData
        (
            mesh,
            patchIDs,
            cellNLayers,
            faceWantedThickness,
            faceRealThickness
        );
    }
}


// ************************************************************************* //
