/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::snappyLayerDriver

Description
    All to do with adding layers

SourceFiles
    snappyLayerDriver.C
    snappyLayerDriverOneByOne.C

\*---------------------------------------------------------------------------*/

#ifndef snappyLayerDriver_H
#define snappyLayerDriver_H

#include "meshRefinement.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class removePoints;
class pointSet;
class motionSmoother;
class addPatchCellLayer;
class faceSet;
class layerParameters;
class externalDisplacementMeshMover;

/*---------------------------------------------------------------------------*\
                           Class snappyLayerDriver Declaration
\*---------------------------------------------------------------------------*/

class snappyLayerDriver
{
public:

    // Public data types

        //- Extrusion controls
        enum extrudeMode
        {
            NOEXTRUDE,      //!< Do not extrude. No layers added.
            EXTRUDE,        //!< Extrude
            EXTRUDEREMOVE   //!< Extrude but afterwards remove added
                            //!< faces locally
        };

private:

    // Private classes

        //- Combine operator class to combine normal with other normal.
        class nomalsCombine
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (y != point::max)
                {
                    if (x == point::max)
                    {
                        x = y;
                    }
                    else
                    {
                        x *= (x&y);
                    }
                }
            }
        };


    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- From surface region to patch
        const labelList globalToMasterPatch_;

        //- From surface region to patch
        const labelList globalToSlavePatch_;

        //- Are we operating in test mode?
        const bool dryRun_;


    // Private Member Functions

        // Layers

            //- For debugging: Dump displacement to .obj files
            static void dumpDisplacement
            (
                const fileName&,
                const indirectPrimitivePatch&,
                const vectorField&,
                const List<extrudeMode>&
            );

            //- Average point wise data to face wise
            static tmp<scalarField> avgPointData
            (
                const indirectPrimitivePatch&,
                const scalarField& pointFld
            );

            //- Check that primitivePatch is not multiply connected.
            //  Collect non-manifold points in pointSet.
            static void checkManifold
            (
                const indirectPrimitivePatch&,
                pointSet& nonManifoldPoints
            );

            //- Check that mesh outside is not multiply connected.
            void checkMeshManifold() const;


            // Static extrusion setup

                //- Unset extrusion on point. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const label patchPointi,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Unset extrusion on face. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const face& localFace,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Truncate index in face
                static label constrainFp(const label sz, const label fp);

                //- Count common points between face and its neighbours
                void countCommonPoints
                (
                    const indirectPrimitivePatch&,
                    const label facei,
                    Map<label>&
                ) const;

                //- Check if any common points form single string. Return
                //  false if not.
                bool checkCommonOrder
                (
                    const label nCommon,
                    const face& curFace,
                    const face& nbFace
                ) const;

                //- Check if any common points form single string; unmark
                //  points on face if not
                void checkCommonOrder
                (
                    const indirectPrimitivePatch& pp,
                    const label facei,
                    const Map<label>& nCommonPoints,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Check if any common points form single string; unmark
                //  points on face if not
                void handleNonStringConnected
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion at non-manifold points.
                void handleNonManifolds
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const labelListList& edgeGlobalFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on feature edges. Assumes non-manifold
                //  edges already handled.
                void handleFeatureAngle
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const scalar minAngle,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on warped faces
                void handleWarpedFaces
                (
                    const indirectPrimitivePatch& pp,
                    const scalar faceRatio,
                    const boolList& relativeSizes,
                    const scalar edge0Len,
                    const labelList& cellLevel,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Determine the number of layers per point from the number of
                //  layers per surface.
                static void setNumLayers
                (
                    meshRefinement& meshRefiner,
                    const labelList& patchToNLayers,
                    const labelList& patchIDs,
                    const indirectPrimitivePatch& pp,

                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus,
                    label& nIdealAddedCells
                );

                //- Determine number of layers per point; include static checks
                //- on invalid extrusion (e.g. non-manifold)
                label setPointNumLayers
                (
                    const layerParameters& layerParams,

                    const labelList& numLayers,
                    const labelList& patchIDs,
                    const indirectPrimitivePatch& pp,
                    const labelListList& edgeGlobalFaces,

                    vectorField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>&
                ) const;
                autoPtr<externalDisplacementMeshMover> makeMeshMover
                (
                    const layerParameters& layerParams,
                    const dictionary& motionDict,
                    const labelList& internalFaceZones,
                    const scalarIOField& minThickness,
                    pointVectorField& displacement
                ) const;
                void updatePatch
                (
                    const labelList& patchIDs,
                    const mapPolyMesh& map,
                    autoPtr<indirectPrimitivePatch>& pp,
                    labelList& newToOldPatchPoints
                ) const;


                //- Helper function to make a pointVectorField with correct
                //  bcs for layer addition:
                //  - numLayers > 0         : fixedValue
                //  - numLayers == 0        : fixedValue (always zero)
                //  - processor             : calculated (so free to move)
                //  - cyclic/wedge/symmetry : slip
                //  - other                 : slip
                static tmp<pointVectorField> makeLayerDisplacementField
                (
                    const pointMesh& pMesh,
                    const labelList& numLayers
                );

                //- Grow no-extrusion layer.
                void growNoExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Re-merge points/faces on faceZones. Opposite of
                //- dupFaceZonePoints above
                void mergeFaceZonePoints
                (
                    const labelList& pointToMaster,
                    labelList& cellNLayers,
                    scalarField& faceRealThickness,
                    scalarField& faceWantedThickness
                );

                //- Calculate pointwise wanted and minimum thickness.
                //  thickness: wanted thickness
                //  minthickness: when to give up and not extrude
                //  Gets per patch parameters and determine pp pointwise
                //  parameters.
                void calculateLayerThickness
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& patchIDs,
                    const layerParameters& layerParams,
                    const labelList& cellLevel,
                    const labelList& patchNLayers,
                    const scalar edge0Len,

                    scalarField& thickness,
                    scalarField& minThickness,
                    scalarField& expansionRatio
                ) const;


            // Extrusion execution

                //- Synchronize displacement among coupled patches.
                static void syncPatchDisplacement
                (
                    const fvMesh& mesh,
                    const indirectPrimitivePatch& pp,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Get nearest point on surface to snap to
                void getPatchDisplacement
                (
                    const indirectPrimitivePatch& pp,
                    const scalarField& thickness,
                    const scalarField& minThickness,
                    const scalarField& expansionRatio,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- For truncateDisplacement: find strings of edges
                bool sameEdgeNeighbour
                (
                    const labelListList& globalEdgeFaces,
                    const label myGlobalFacei,
                    const label nbrGlobFacei,
                    const label edgeI
                ) const;

                //- For truncateDisplacement: find strings of edges
                void getVertexString
                (
                    const indirectPrimitivePatch& pp,
                    const labelListList& globalEdgeFaces,
                    const label facei,
                    const label edgeI,
                    const label myGlobFacei,
                    const label nbrGlobFacei,
                    DynamicList<label>& vertices
                ) const;

                //- Truncates displacement
                // - for all patchFaces in the faceset displacement gets set
                //   to zero
                // - all displacement < minThickness gets set to zero
                // - all non-consecutive extrusions get set to 0
                label truncateDisplacement
                (
                    const globalIndex& globalFaces,
                    const labelListList& edgeGlobalFaces,
                    const indirectPrimitivePatch& pp,
                    const scalarField& minThickness,
                    const faceSet& illegalPatchFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Setup layer information (at points and faces) to
                //  modify mesh topology in
                //  regions where layer mesh terminates. Guarantees an
                //  optional slow decreasing of the number of layers.
                //  Returns the number of layers per face and per point
                //  to go into the actual layer addition engine.
                void setupLayerInfoTruncation
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& patchNLayers,
                    const List<extrudeMode>& extrudeStatus,
                    const label nBufferCellsNoExtrude,
                    labelList& nPatchPointLayers,
                    labelList& nPatchFaceLayers
                ) const;

                //- Does any of the cells use a face from faces?
                static bool cellsUseFace
                (
                    const polyMesh& mesh,
                    const labelList& cellLabels,
                    const labelHashSet& faces
                );

                //- Checks the newly added cells and locally unmarks points
                //  so they will not get extruded next time round. Returns
                //  global number of unmarked points (0 if all was fine)
                static label checkAndUnmark
                (
                    const addPatchCellLayer& addLayer,
                    const dictionary& motionDict,
                    const bool additionalReporting,
                    const List<labelPair>& baffles,
                    const indirectPrimitivePatch& pp,
                    const fvMesh&,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Count global number of extruded faces
                static label countExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    const List<extrudeMode>& extrudeStatus
                );

                //- After adding to mesh get the new baffles
                static List<labelPair> getBafflesOnAddedMesh
                (
                    const polyMesh& mesh,
                    const labelList& newToOldFaces,
                    const List<labelPair>& baffles
                );

                //- Collect layer faces and layer cells into bools
                //  for ease of handling
                static void getLayerCellsFaces
                (
                    const polyMesh&,
                    const addPatchCellLayer&,
                    const scalarField& oldRealThickness,

                    labelList& cellStatus,
                    scalarField& faceRealThickness
                );

                //- Print layer coverage table
                void printLayerData
                (
                    const fvMesh& mesh,
                    const labelList& patchIDs,
                    const labelList& cellNLayers,
                    const scalarField& faceWantedThickness,
                    const scalarField& faceRealThickness,
                    const layerParameters& layerParams
                ) const;

                //- Write cellSet,faceSet for layers
                bool writeLayerSets
                (
                    const fvMesh& mesh,
                    const labelList& cellNLayers,
                    const scalarField& faceRealThickness
                ) const;

                //- Write volFields,cellSet,faceSet for layers depending
                //  on write level
                bool writeLayerData
                (
                    const fvMesh& mesh,
                    const labelList& patchIDs,
                    const labelList& cellNLayers,
                    const scalarField& faceWantedThickness,
                    const scalarField& faceRealThickness
                ) const;


            // Mesh shrinking (to create space for layers)

                //- Average field (over all subset of mesh points) by
                //  summing contribution from edges. Global parallel since only
                //  does master edges for coupled edges.
                template<class Type>
                static void averageNeighbours
                (
                    const polyMesh& mesh,
                    const bitSet& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    const scalarField& invSumWeight,
                    const Field<Type>& data,
                    Field<Type>& average
                );

                //- Calculate inverse sum of edge weights (currently always 1.0)
                void sumWeights
                (
                    const bitSet& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    scalarField& invSumWeight
                ) const;

                //- Smooth scalar field on patch
                void smoothField
                (
                    const motionSmoother& meshMover,
                    const bitSet& isMasterPoint,
                    const bitSet& isMasterEdge,
                    const labelList& meshEdges,
                    const scalarField& fieldMin,
                    const label nSmoothDisp,
                    scalarField& field
                ) const;

                //- Smooth normals on patch.
                void smoothPatchNormals
                (
                    const motionSmoother& meshMover,
                    const bitSet& isMasterPoint,
                    const bitSet& isMasterEdge,
                    const labelList& meshEdges,
                    const label nSmoothDisp,
                    pointField& normals
                ) const;

                //- Smooth normals in interior.
                void smoothNormals
                (
                    const label nSmoothDisp,
                    const bitSet& isMasterPoint,
                    const bitSet& isMasterEdge,
                    const labelList& fixedPoints,
                    pointVectorField& normals
                ) const;

                //- Stop layer growth where mesh wraps around edge with a
                //  large feature angle
                void handleFeatureAngleLayerTerminations
                (
                    const scalar minCos,
                    const bitSet& isMasterPoint,
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,

                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    label& nPointCounter
                ) const;

                //- Find isolated islands (points, edges and faces and
                // layer terminations)
                // in the layer mesh and stop any layer growth at these points.
                void findIsolatedRegions
                (
                    const scalar minCosLayerTermination,
                    const bitSet& isMasterPoint,
                    const bitSet& isMasterEdge,
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const scalarField& minThickness,

                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers
                ) const;


        //- No copy construct
        snappyLayerDriver(const snappyLayerDriver&) = delete;

        //- No copy assignment
        void operator=(const snappyLayerDriver&) = delete;


public:

    //- Runtime type information
    ClassName("snappyLayerDriver");

    // Constructors

        //- Construct from components
        snappyLayerDriver
        (
            meshRefinement& meshRefiner,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch,
            const bool dryRun = false
        );


    // Member Functions

            //- Merge patch faces on same cell.
            void mergePatchFacesUndo
            (
                const layerParameters& layerParams,
                const dictionary& motionDict,
                const meshRefinement::FaceMergeType mergeType
            );

            void addLayers
            (
                const layerParameters& layerParams,
                const label nLayerIter,

                const dictionary& motionDict,
                const label nRelaxedIter,
                const label nAllowableErrors,

                const labelList& patchIDs,
                const labelList& internalFaceZones,
                const List<labelPair>& baffles,
                const labelList& numLayers,
                const label nIdealTotAddedCells,

                const globalIndex& globalFaces,
                indirectPrimitivePatch& pp,
                const labelListList& edgeGlobalFaces,
                const labelList& edgePatchID,
                const labelList& edgeZoneID,
                const boolList& edgeFlip,
                const labelList& inflateFaceID,
                const scalarField& thickness,
                const scalarIOField& minThickness,
                const scalarField& expansionRatio,
                vectorField& patchDisp,
                labelList& patchNLayers,
                List<extrudeMode>& extrudeStatus,
                polyTopoChange& savedMeshMod,
                labelList& cellNLayers,
                scalarField& faceRealThickness
            );

            //- Add cell layers
            void addLayers
            (
                const layerParameters& layerParams,
                const dictionary& motionDict,
                const labelList& patchIDs,
                const label nAllowableErrors,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- For debugging. Can be removed.
            void addLayersSinglePass
            (
                const layerParameters& layerParams,
                const dictionary& motionDict,
                const labelList& patchIDs,
                const label nAllowableErrors,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Add layers according to the dictionary settings
            void doLayers
            (
                const dictionary& shrinkDict,
                const dictionary& motionDict,
                const layerParameters& layerParams,
                const meshRefinement::FaceMergeType mergeType,
                const bool preBalance,              // balance before adding?
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Helper: see what zones and patches edges should be extruded into
            static void determineSidePatches
            (
                meshRefinement& meshRefiner,
                const globalIndex& globalFaces,
                const labelListList& edgeGlobalFaces,
                const indirectPrimitivePatch& pp,

                labelList& edgePatchID,
                labelList& edgeZoneID,
                boolList& edgeFlip,
                labelList& inflateFaceID
            );

            //- Duplicate points on faceZones with layers. Re-used when adding
            //- buffer layers. Can be made private again once multi-side
            //- layer addition working.
            static autoPtr<mapPolyMesh> dupFaceZonePoints
            (
                meshRefinement& meshRefiner,
                const labelList& patchIDs,  // patch indices
                const labelList& numLayers, // number of layers per patch
                List<labelPair> baffles,
                labelList& pointToMaster
            );

            //- Map numbering after adding cell layers
            static void mapFaceZonePoints
            (
                meshRefinement& meshRefiner,
                const mapPolyMesh& map,
                labelPairList& baffles,
                labelList& pointToMaster
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "snappyLayerDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
