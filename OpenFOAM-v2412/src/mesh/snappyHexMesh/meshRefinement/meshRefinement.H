/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::meshRefinement

Description
    Helper class which maintains intersections of (changing) mesh with
    (static) surfaces.

    Maintains
    - per face any intersections of the cc-cc segment with any of the surfaces

SourceFiles
    meshRefinement.C
    meshRefinementBaffles.C
    meshRefinementGapRefine.C
    meshRefinementMerge.C
    meshRefinementProblemCells.C
    meshRefinementRefine.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_meshRefinement_H
#define Foam_meshRefinement_H

#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "autoPtr.H"
#include "labelPairHashes.H"
#include "indirectPrimitivePatch.H"
#include "pointFieldsFwd.H"
#include "Tuple2.H"
#include "pointIndexHit.H"
#include "wordPairHashes.H"
#include "surfaceZonesInfo.H"
#include "volumeType.H"
#include "DynamicField.H"
#include "coordSetWriter.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class mapDistributePolyMesh;
class mapDistribute;
class decompositionMethod;
class refinementSurfaces;
class refinementFeatures;
class shellSurfaces;
class removeCells;
class fvMeshDistribute;
class removePoints;
class localPointRegion;
class snapParameters;

/*---------------------------------------------------------------------------*\
                       Class meshRefinement Declaration
\*---------------------------------------------------------------------------*/

class meshRefinement
{
public:

    // Public Data Types

        //- Enumeration for how to operate
        enum MeshType
        {
            CASTELLATED,    // split-hex refinement, snapping, layer addition
            CASTELLATEDBUFFERLAYER, // same but add buffer layers before
                                    // snapping
            CASTELLATEDBUFFERLAYER2 // experimental. wip.
        };

        static const Enum<MeshType> MeshTypeNames;


        //- Enumeration for what to debug. Used as a bit-pattern.
        enum debugType
        {
            MESH = (1 << 0),
            OBJINTERSECTIONS = (1 << 1),
            FEATURESEEDS = (1 << 2),
            ATTRACTION = (1 << 3),
            LAYERINFO = (1 << 4)
        };

        static const Enum<debugType> debugTypeNames;

        ////- Enumeration for what to output. Used as a bit-pattern.
        //enum outputType
        //{
        //    OUTPUTLAYERINFO = (1 << 0)
        //};
        //
        //static const Enum<outputType> outputTypeNames;

        //- Enumeration for what to write. Used as a bit-pattern.
        enum writeType
        {
            WRITEMESH = (1 << 0),
            NOWRITEREFINEMENT = (1 << 1),
            WRITELEVELS = (1 << 2),
            WRITELAYERSETS = (1 << 3),
            WRITELAYERFIELDS = (1 << 4)
        };

        static const Enum<writeType> writeTypeNames;

        //- Enumeration for how userdata is to be mapped upon refinement.
        enum mapType
        {
            MASTERONLY = 1, //!< maintain master only
            KEEPALL = 2,    //!< have slaves (upon refinement) from master
            REMOVE = 4      //!< set value to -1 any face that was refined
        };

        //- Enumeration for what to do with co-planar patch faces on a single
        //  cell
        enum FaceMergeType
        {
            NONE,           // no merging
            GEOMETRIC,      // use feature angle
            IGNOREPATCH     // use feature angle, allow merging of different
                            // patches
        };


private:

    // Static Data Members

        //- Control of writing level
        static writeType writeLevel_;

        ////- Control of output/log level
        //static outputType outputLevel_;


    // Private Data

        //- Reference to mesh
        fvMesh& mesh_;

        //- Tolerance used for sorting coordinates (used in 'less' routine)
        const scalar mergeDistance_;

        //- Overwrite the mesh?
        const bool overwrite_;

        //- Instance of mesh upon construction. Used when in overwrite_ mode.
        const word oldInstance_;

        //- All surface-intersection interaction
        const refinementSurfaces& surfaces_;

        //- All feature-edge interaction
        const refinementFeatures& features_;

        //- All shell-refinement interaction
        const shellSurfaces& shells_;

        //- All limit-refinement interaction
        const shellSurfaces& limitShells_;

        //- How to generate mesh
        const MeshType meshType_;

        //- Are we operating in test mode?
        const bool dryRun_;


        //- Refinement engine
        hexRef8 meshCutter_;

        //- Per cc-cc vector the index of the surface hit
        labelIOList surfaceIndex_;


        // For baffle merging

            //- Original patch for baffle faces that used to be on
            //  coupled patches
            Map<label> faceToCoupledPatch_;


        //- User supplied face based data.
        List<Tuple2<mapType, labelList>> userFaceData_;

        //- Meshed patches - are treated differently. Stored as wordList since
        //  order changes.
        wordList meshedPatches_;

        //- FaceZone to master patch name
        HashTable<word> faceZoneToMasterPatch_;

        //- FaceZone to slave patch name
        HashTable<word> faceZoneToSlavePatch_;

        //- FaceZone to method to handle faces
        HashTable<surfaceZonesInfo::faceZoneType> faceZoneToType_;


    // Private Member Functions

        //- Add patchfield of given type to all fields on mesh
        template<class GeoField>
        static void addPatchFields(fvMesh&, const word& patchFieldType);

        //- Reorder patchfields of all fields on mesh
        template<class GeoField>
        static void reorderPatchFields(fvMesh&, const labelList& oldToNew);

        //- Find out which faces have changed given cells (old mesh labels)
        //  that were marked for refinement.
        static labelList getChangedFaces
        (
            const mapPolyMesh&,
            const labelList& oldCellsToRefine
        );

        // Debug: count number of master-faces (and do some checking
        // for consistency)
        label globalFaceCount(const labelList& elems) const;

        //- Calculate coupled boundary end vector and refinement level
        void calcNeighbourData(labelList& neiLevel, pointField& neiCc) const;

        //- Calculate rays from cell-centre to cell-centre and corresponding
        //  min cell refinement level
        void calcCellCellRays
        (
            const pointField& neiCc,
            const labelList& neiLevel,
            const labelList& testFaces,
            pointField& start,
            pointField& end,
            labelList& minLevel
        ) const;

        //- Get cells which are inside any closed surface. Note that
        //  all closed surfaces
        //  will have already been oriented to have keepPoint outside.
        labelList getInsideCells(const word&) const;

        //- Do all to remove inside cells
        autoPtr<mapPolyMesh> removeInsideCells
        (
            const string& msg,
            const label exposedPatchi
        );


        // Refinement candidate selection

            //- Mark cell for refinement (if not already marked). Return false
            // if refinelimit hit. Keeps running count (in nRefine) of cells
            // marked for refinement
            static bool markForRefine
            (
                const label markValue,
                const label nAllowRefine,
                label& cellValue,
                label& nRefine
            );

            //- Mark every cell with level of feature passing through it
            //  (or -1 if not passed through). Uses tracking.
            void markFeatureCellLevel
            (
                const pointField& keepPoints,
                labelList& maxFeatureLevel
            ) const;

            //- Calculate list of cells to refine based on intersection of
            //  features.
            label markFeatureRefinement
            (
                const pointField& keepPoints,
                const label nAllowRefine,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for distance-to-feature based refinement.
            label markInternalDistanceToFeatureRefinement
            (
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for refinement-shells based refinement.
            label markInternalRefinement
            (
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Unmark cells for refinement based on limit-shells. Return number
            //  of limited cells.
            label unmarkInternalRefinement
            (
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Collect faces that are intersected and whose neighbours aren't
            //  yet marked for refinement.
            labelList getRefineCandidateFaces
            (
                const labelList& refineCell
            ) const;

            //- Mark cells for surface intersection based refinement.
            label markSurfaceRefinement
            (
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Collect cells intersected by the surface that are candidates
            //  for gap checking. Used inside markSurfaceGapRefinement
            void collectGapCandidates
            (
                const shellSurfaces& shells,
                const labelList& testFaces,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& cellToCompact,
                labelList& bFaceToCompact,
                List<FixedList<label, 3>>& shellGapInfo,
                List<volumeType>& shellGapMode
            ) const;
            void collectGapCells
            (
                const scalar planarCos,

                const List<FixedList<label, 3>>& extendedGapLevel,
                const List<volumeType>& extendedGapMode,
                const labelList& testFaces,
                const pointField& start,
                const pointField& end,

                const labelList& cellToCompact,
                const labelList& bFaceToCompact,
                const List<FixedList<label, 3>>& shellGapInfo,
                const List<volumeType>& shellGapMode,

                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;


            //- Mark cells intersected by the surface if they are inside
            //  close gaps
            label markSurfaceGapRefinement
            (
                const scalar planarCos,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Generate single ray from nearPoint in direction of nearNormal
            label generateRays
            (
                const point& nearPoint,
                const vector& nearNormal,
                const FixedList<label, 3>& gapInfo,
                const volumeType& mode,

                const label cLevel,

                DynamicField<point>& start,
                DynamicField<point>& end
            ) const;

            //- Generate pairs of rays through cell centre
            //  Each ray pair has start, end, and expected gap size
            label generateRays
            (
                const bool useSurfaceNormal,

                const point& nearPoint,
                const vector& nearNormal,
                const FixedList<label, 3>& gapInfo,
                const volumeType& mode,

                const point& cc,
                const label cLevel,

                DynamicField<point>& start,
                DynamicField<point>& end,
                DynamicField<scalar>& gapSize,

                DynamicField<point>& start2,
                DynamicField<point>& end2,
                DynamicField<scalar>& gapSize2
            ) const;

            //- Select candidate cells (cells inside a shell with gapLevel
            //  specified)
            void selectGapCandidates
            (
                const labelList& refineCell,
                const label nRefine,

                labelList& cellMap,
                labelList& gapShell,
                List<FixedList<label, 3>>& shellGapInfo,
                List<volumeType>& shellGapMode
            ) const;

            //- Merge gap information coming from shell and from surface
            //  (surface wins)
            void mergeGapInfo
            (
                const FixedList<label, 3>& shellGapInfo,
                const volumeType shellGapMode,
                const FixedList<label, 3>& surfGapInfo,
                const volumeType surfGapMode,
                FixedList<label, 3>& gapInfo,
                volumeType& gapMode
            ) const;

            //- Mark cells for non-surface intersection based gap refinement
            label markInternalGapRefinement
            (
                const scalar planarCos,
                const bool spreadGapSize,
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine,
                labelList& numGapCells,
                scalarField& gapSize
            ) const;

            //- Refine cells containing small gaps
            label markSmallFeatureRefinement
            (
                const scalar planarCos,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Refine cells containing triangles with high refinement level
            //  (currently due to high curvature or being inside shell with
            //   high level)
            label markSurfaceFieldRefinement
            (
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Helper: count number of normals1 that are in normals2
            label countMatches
            (
                const List<point>& normals1,
                const List<point>& normals2,
                const scalar tol = 1e-6
            ) const;

            //- Detect if two intersection points are high curvature (w.r.t.
            //  lengthscale
            bool highCurvature
            (
                const scalar minCosAngle,
                const scalar lengthScale,
                const point& p0,
                const vector& n0,
                const point& p1,
                const vector& n1
            ) const;

            //- Mark cells for surface curvature based refinement. Marks if
            //  local curvature > curvature or if on different regions
            //  (markDifferingRegions)
            label markSurfaceCurvatureRefinement
            (
                const scalar curvature,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cell if local a gap topology or
            bool checkProximity
            (
                const scalar planarCos,
                const label nAllowRefine,

                const label surfaceLevel,
                const vector& surfaceLocation,
                const vector& surfaceNormal,

                const label celli,

                label& cellMaxLevel,
                vector& cellMaxLocation,
                vector& cellMaxNormal,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for surface proximity based refinement.
            label markProximityRefinement
            (
                const scalar curvature,

                // Per region the min and max cell level
                const labelList& surfaceMinLevel,
                const labelList& surfaceMaxLevel,

                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark faces for additional proximity based testing. Extends
            //  testFaces
            void markProximityCandidateFaces
            (
                const labelList& blockedSurfaces,
                const List<scalarList>& regionToBlockSize,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& testFaces
            ) const;

            //- Mark cells for surface proximity based refinement. Uses
            //  ray-shooting from markInternalGapRefinement.
            label markFakeGapRefinement
            (
                const scalar planarCos,
                const labelList& blockedSurfaces,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;


        // Baffle handling

            //- Get faces to repatch. Returns map from face to patch.
            Map<labelPair> getZoneBafflePatches
            (
                const bool allowBoundary,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch
            ) const;

            //- Calculate intersections. Return per face -1 or the global
            //  surface region
            void getIntersections
            (
                const labelList& surfacesToTest,
                const pointField& neiCc,
                const labelList& testFaces,

                labelList& globalRegion1,
                labelList& globalRegion2
            ) const;

            //- Calculate intersections on zoned faces. Return per face -1
            //  or the global region of the surface and the orientation
            //  w.r.t. surface
            void getIntersections
            (
                const labelList& surfacesToTest,
                const pointField& neiCc,
                const labelList& testFaces,

                labelList& namedSurfaceRegion,
                bitSet& posOrientation
            ) const;

            //- Determine patches for baffles
            void getBafflePatches
            (
                const label nErodeCellZones,
                const labelList& globalToMasterPatch,
                const pointField& locationsInMesh,
                const wordList& regionsInMesh,
                const pointField& locationsOutsideMesh,
                const bool exitIfLeakPath,
                const refPtr<coordSetWriter>& leakPathFormatter,
                refPtr<surfaceWriter>& surfFormatter,

                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& ownPatch,
                labelList& neiPatch
            ) const;

            autoPtr<mapPolyMesh> splitMesh
            (
                const label nBufferLayers,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                labelList& cellRegion,
                labelList& ownPatch,
                labelList& neiPatch
            );

            //- Repatches external face or creates baffle for internal face
            //  with user specified patches (might be different for both sides).
            //  Returns label of added face.
            label createBaffle
            (
                const label facei,
                const label ownPatch,
                const label neiPatch,
                polyTopoChange& meshMod
            ) const;

            //- Write leak path
            static fileName writeLeakPath
            (
                const polyMesh& mesh,
                const pointField& locationsInMesh,
                const pointField& locationsOutsideMesh,
                const boolList& blockedFace,
                coordSetWriter& leakPathWriter
            );


        // Problem cell handling

            //- Helper function to mark face as being on 'boundary'. Used by
            //  markFacesOnProblemCells
            void markBoundaryFace
            (
                const label facei,
                boolList& isBoundaryFace,
                boolList& isBoundaryEdge,
                boolList& isBoundaryPoint
            ) const;

            void findNearest
            (
                const labelList& meshFaces,
                List<pointIndexHit>& nearestInfo,
                labelList& nearestSurface,
                labelList& nearestRegion,
                vectorField& nearestNormal
            ) const;

            Map<label> findEdgeConnectedProblemCells
            (
                const scalarField& perpendicularAngle,
                const labelList&
            ) const;

            bool isCollapsedFace
            (
                const pointField&,
                const pointField& neiCc,
                const scalar minFaceArea,
                const scalar maxNonOrtho,
                const label facei
            ) const;

            bool isCollapsedCell
            (
                const pointField&,
                const scalar volFraction,
                const label celli
            ) const;

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into. If removeEdgeConnectedCells is set
            //  removes cells based on perpendicularAngle.
            void markFacesOnProblemCells
            (
                const dictionary& motionDict,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const labelList& globalToMasterPatch,

                labelList& facePatch,
                labelList& faceZone
            ) const;

            //- Calculates for every face the nearest 'start' face. Any
            //  unreached face does not get set (faceToStart[facei] = -1)
            void nearestFace
            (
                const labelUList& startFaces,
                const bitSet& isBlockedFace,

                autoPtr<mapDistribute>& mapPtr,
                labelList& faceToStart,
                const label nIter = labelMax
            ) const;

            //- Calculates for every face the label of the nearest
            //  patch its zone. Any unreached face (disconnected mesh?) becomes
            //  adaptPatchIDs[0]
            void nearestPatch
            (
                const labelList& adaptPatchIDs,
                labelList& nearestPatch,
                labelList& nearestZone
            ) const;

            //- Returns list with for every face the label of the nearest
            //  patch. Any unreached face (disconnected mesh?) becomes
            //  adaptPatchIDs[0]
            labelList nearestPatch(const labelList& adaptPatchIDs) const;

            //- Returns list with for every face the label of the nearest
            //  (global) region. Any unreached face (disconnected mesh?) becomes
            //  defaultRegion
            labelList nearestIntersection
            (
                const labelList& surfacesToTest,
                const label defaultRegion
            ) const;

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into.
            void markFacesOnProblemCellsGeometric
            (
                const snapParameters& snapParams,
                const dictionary& motionDict,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,

                labelList& facePatch,
                labelList& faceZone
            ) const;


        // Baffle merging

            //- Extract those baffles (duplicate) faces that are on the edge
            //  of a baffle region. These are candidates for merging.
            //  samePatch : include only if both sides on the same patch.
            List<labelPair> freeStandingBaffles
            (
                const List<labelPair>&,
                const scalar freeStandingAngle,
                const bool samePatch
            ) const;


        // Zone handling

            //- Finds zone per cell for cells inside closed named surfaces.
            //  (uses geometric test for insideness)
            //  Adapts namedSurfaceRegion so all faces on boundary of cellZone
            //  have corresponding faceZone.
            void findCellZoneGeometric
            (
                const pointField& neiCc,
                const labelList& closedNamedSurfaces,
                labelList& namedSurfaceRegion,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Finds zone per cell for cells inside region for which name
            //  is specified.
            void findCellZoneInsideWalk
            (
                const pointField& locationsInMesh,
                const labelList& zonesInMesh,
                const labelList& blockedFace, // per face -1 or some index >= 0
                labelList& cellToZone
            ) const;

            //- Finds zone per cell for cells inside region for which name
            //  is specified.
            void findCellZoneInsideWalk
            (
                const pointField& locationsInMesh,
                const wordList& zoneNamesInMesh,
                const labelList& faceToZone, // per face -1 or some index >= 0
                labelList& cellToZone
            ) const;

            //- Determines cell zone from cell region information.
            bool calcRegionToZone
            (
                const label backgroundZoneID,
                const label surfZoneI,
                const label ownRegion,
                const label neiRegion,

                labelList& regionToCellZone
            ) const;

            //- Finds zone per cell. Uses topological walk with all faces
            //  marked in unnamedSurfaceRegion (intersections with unnamed
            //  surfaces) and namedSurfaceRegion (intersections with named
            //  surfaces) regarded as blocked.
            void findCellZoneTopo
            (
                const label backgroundZoneID,
                const pointField& locationsInMesh,
                const labelList& unnamedSurfaceRegion,
                const labelList& namedSurfaceRegion,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Opposite of findCellTopo: finds assigned cell connected to
            //  an unassigned one and puts it in the background zone.
            void erodeCellZone
            (
                const label nErodeCellZones,
                const label backgroundZoneID,
                const labelList& unnamedSurfaceRegion,
                const labelList& namedSurfaceRegion,
                labelList& cellToZone
            ) const;

            //- Variation of findCellZoneTopo: walks from cellZones but ignores
            //  face intersections (unnamedSurfaceRegion). WIP
            void growCellZone
            (
                const label nGrowCellZones,
                const label backgroundZoneID,
                labelList& unnamedSurfaceRegion1,
                labelList& unnamedSurfaceRegion2,
                labelList& namedSurfaceRegion,
                labelList& cellToZone
            ) const;

            //- Make namedSurfaceRegion consistent with cellToZone
            //  - clear out any blocked faces inbetween same cell zone.
            void makeConsistentFaceIndex
            (
                const labelList& zoneToNamedSurface,
                const labelList& cellToZone,
                labelList& namedSurfaceRegion
            ) const;

            //- Calculate cellZone allocation
            void zonify
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
            ) const;

            //- Put cells into cellZone, faces into faceZone
            void zonify
            (
                const bitSet& isMasterFace,
                const labelList& cellToZone,
                const labelList& neiCellZone,
                const labelList& faceToZone,
                const bitSet& meshFlipMap,
                polyTopoChange& meshMod
            ) const;

            //- Allocate faceZoneName
            void allocateInterRegionFaceZone
            (
                const label ownZone,
                const label neiZone,
                wordPairHashTable& zonesToFaceZone,
                LabelPairMap<word>& zoneIDsToFaceZone
            ) const;

            //- Remove any loose standing cells
            void handleSnapProblems
            (
                const snapParameters& snapParams,
                const bool useTopologicalSnapDetection,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch
            );


            // Some patch utilities

            //- Get all faces in faceToZone that have no cellZone on
            //  either side.
            labelList freeStandingBaffleFaces
            (
                const labelList& faceToZone,
                const labelList& cellToZone,
                const labelList& neiCellZone
            ) const;

            //- Determine per patch edge the number of master faces. Used
            //  to detect non-manifold situations.
            void calcPatchNumMasterFaces
            (
                const bitSet& isMasterFace,
                const indirectPrimitivePatch& patch,
                labelList& nMasterFaces
            ) const;

            //- Determine per patch face the (singly-) connected zone it
            //  is in. Return overall number of zones.
            label markPatchZones
            (
                const indirectPrimitivePatch& patch,
                const labelList& nMasterFaces,
                labelList& faceToZone
            ) const;

            //- Make faces consistent.
            void consistentOrientation
            (
                const bitSet& isMasterFace,
                const indirectPrimitivePatch& patch,
                const labelList& nMasterFaces,
                const labelList& faceToZone,
                const Map<label>& zoneToOrientation,
                bitSet& meshFlipMap
            ) const;


        //- No copy construct
        meshRefinement(const meshRefinement&) = delete;

        //- No copy assignment
        void operator=(const meshRefinement&) = delete;

public:

    //- Runtime type information
    ClassName("meshRefinement");


    // Constructors

        //- Construct from components
        meshRefinement
        (
            fvMesh& mesh,
            const scalar mergeDistance,
            const bool overwrite,
            const refinementSurfaces&,
            const refinementFeatures&,
            const shellSurfaces&,   // omnidirectional refinement
            const shellSurfaces&,   // limit refinement
            const labelUList& checkFaces,   // initial faces to check
            const MeshType meshType,
            const bool dryRun
        );


    // Member Functions

        // Access

            //- Reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }
            fvMesh& mesh()
            {
                return mesh_;
            }

            scalar mergeDistance() const
            {
                return mergeDistance_;
            }

            //- Overwrite the mesh?
            bool overwrite() const
            {
                return overwrite_;
            }

            //- (points)instance of mesh upon construction
            const word& oldInstance() const
            {
                return oldInstance_;
            }

            //- Reference to surface search engines
            const refinementSurfaces& surfaces() const
            {
                return surfaces_;
            }

            //- Reference to feature edge mesh
            const refinementFeatures& features() const
            {
                return features_;
            }

            //- Reference to refinement shells (regions)
            const shellSurfaces& shells() const
            {
                return shells_;
            }

            //- Reference to limit shells (regions)
            const shellSurfaces& limitShells() const
            {
                return limitShells_;
            }

            //- Reference to meshcutting engine
            const hexRef8& meshCutter() const
            {
                return meshCutter_;
            }

            //- Mode of meshing
            MeshType meshType() const
            {
                return meshType_;
            }

            //- Per start-end edge the index of the surface hit
            const labelList& surfaceIndex() const;

            labelList& surfaceIndex();

            //- For faces originating from processor faces store the original
            //  patch
            const Map<label>& faceToCoupledPatch() const
            {
                return faceToCoupledPatch_;
            }

            //- Additional face data that is maintained across
            //  topo changes. Every entry is a list over all faces.
            //  Bit of a hack. Additional flag to say whether to maintain master
            //  only (false) or increase set to account for face-from-face.
            const List<Tuple2<mapType, labelList>>& userFaceData() const
            {
                return userFaceData_;
            }

            List<Tuple2<mapType, labelList>>& userFaceData()
            {
                return userFaceData_;
            }


        // Other

            //- Count number of intersections (local)
            label countHits() const;

            //- Redecompose according to cell count
            //  keepZoneFaces : find all faceZones from zoned surfaces and keep
            //                  owner and neighbour together
            //  keepBaffles   : find all baffles and keep them together
            //  singleProcPoints : all faces using point are kept together
            autoPtr<mapDistributePolyMesh> balance
            (
                const bool keepZoneFaces,
                const bool keepBaffles,
                const labelList& singleProcPoints,
                const scalarField& cellWeights,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Get faces with intersection.
            labelList intersectedFaces() const;

            //- Get points on surfaces with intersection and boundary faces.
            labelList intersectedPoints() const;

            //- Create patch from set of patches
            static autoPtr<indirectPrimitivePatch> makePatch
            (
                const polyMesh&,
                const labelList&
            );

            //- Helper function to make a pointVectorField with correct
            //  bcs for mesh movement:
            //  - adaptPatchIDs         : fixedValue
            //  - processor             : calculated (so free to move)
            //  - cyclic/wedge/symmetry : slip
            //  - other                 : slip
            static tmp<pointVectorField> makeDisplacementField
            (
                const pointMesh& pMesh,
                const labelList& adaptPatchIDs
            );

            //- Helper function: check that face zones are synced
            static void checkCoupledFaceZones(const polyMesh&);

            //- Helper: calculate edge weights (1/length)
            static void calculateEdgeWeights
            (
                const polyMesh& mesh,
                const bitSet& isMasterEdge,
                const labelList& meshPoints,
                const edgeList& edges,
                scalarField& edgeWeights,
                scalarField& invSumWeight
            );

            //- Helper: weighted sum (over all subset of mesh points) by
            //  summing contribution from (master) edges
            template<class Type>
            static void weightedSum
            (
                const polyMesh& mesh,
                const bitSet& isMasterEdge,
                const labelList& meshPoints,
                const edgeList& edges,
                const scalarField& edgeWeights,
                const Field<Type>& data,
                Field<Type>& sum
            );


        // Refinement

            //- Is local topology a small gap?
            bool isGap
            (
                const scalar,
                const vector&,
                const vector&,
                const vector&,
                const vector&
            ) const;

            //- Is local topology a small gap normal to the test vector
            bool isNormalGap
            (
                const scalar planarCos,
                const label level0,
                const vector& point0,
                const vector& normal0,
                const label level1,
                const vector& point1,
                const vector& normal1
            ) const;

            //- Calculate list of cells to refine.
            labelList refineCandidates
            (
                const pointField& keepPoints,
                const scalar curvature,
                const scalar planarAngle,

                const bool featureRefinement,
                const bool featureDistanceRefinement,
                const bool internalRefinement,
                const bool surfaceRefinement,
                const bool curvatureRefinement,
                const bool smallFeatureRefinement,
                const bool gapRefinement,
                const bool bigGapRefinement,
                const bool spreadGapSize,
                const label maxGlobalCells,
                const label maxLocalCells
            ) const;


            // Blocking cells

                //- Mark faces on interface between set and rest
                //  (and same cell level)
                void markOutsideFaces
                (
                    const labelList& cellLevel,
                    const labelList& neiLevel,
                    const labelList& refineCell,
                    bitSet& isOutsideFace
                ) const;

                //- Count number of faces on cell that are in set
                label countFaceDirs
                (
                    const bitSet& isOutsideFace,
                    const label celli
                ) const;

                //- Add one layer of cells to set
                void growSet
                (
                    const labelList& neiLevel,
                    const bitSet& isOutsideFace,
                    labelList& refineCell,
                    label& nRefine
                ) const;

                //- Detect gapRefinement cells and remove them
                autoPtr<mapPolyMesh> removeGapCells
                (
                    const scalar planarAngle,
                    const labelList& minSurfaceLevel,
                    const labelList& globalToMasterPatch,
                    const label growIter
                );

            // Blocking leaks (by blocking cells)

                //- Faces currently on boundary or intersected by surface
                void selectIntersectedFaces
                (
                    const labelList& surfaces,
                    boolList& isBlockedFace
                ) const;

                //- Return list of cells to block by walking from the seedCells
                //  until reaching a leak face
                labelList detectLeakCells
                (
                    const boolList& isBlockedFace,
                    const labelList& leakFaces,
                    const labelList& seedCells
                ) const;

                //- Remove minimum amount of cells to break any leak from
                //  inside to outside
                autoPtr<mapPolyMesh> removeLeakCells
                (
                    const labelList& globalToMasterPatch,
                    const labelList& globalToSlavePatch,
                    const pointField& locationsInMesh,
                    const pointField& locationsOutsideMesh,
                    const labelList& selectedSurfaces
                );

                //- Baffle faces to break any leak from inside to outside
                autoPtr<mapPolyMesh> blockLeakFaces
                (
                    const labelList& globalToMasterPatch,
                    const labelList& globalToSlavePatch,
                    const pointField& locationsInMesh,
                    const wordList& zonesInMesh,
                    const pointField& locationsOutsideMesh,
                    const labelList& selectedSurfaces
                );


            //- Refine some cells
            autoPtr<mapPolyMesh> refine(const labelList& cellsToRefine);

            //- Balance the mesh
            autoPtr<mapDistributePolyMesh> balance
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                labelList& cellsToRefine,
                const scalar maxLoadUnbalance,
                const label maxCellUnbalance
            );

            //- Refine some cells and rebalance
            autoPtr<mapDistributePolyMesh> refineAndBalance
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const scalar maxLoadUnbalance,
                const label maxCellUnbalance
            );

            //- Balance before refining some cells
            autoPtr<mapDistributePolyMesh> balanceAndRefine
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const scalar maxLoadUnbalance,
                const label maxCellUnbalance
            );

            //- Calculate list of cells to directionally refine
            labelList directionalRefineCandidates
            (
                const label maxGlobalCells,
                const label maxLocalCells,
                const labelList& currentLevel,
                const direction dir
            ) const;

            //- Directionally refine in direction cmpt
            autoPtr<mapPolyMesh> directionalRefine
            (
                const string& msg,
                const direction cmpt,
                const labelList& cellsToRefine
            );


        // Baffle handling

            //- Split off unreachable areas of mesh.
            void baffleAndSplitMesh
            (
                const bool handleSnapProblems,

                // How to remove problem snaps
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
                const wordList& regionsInMesh,
                const pointField& locationsOutsideMesh,
                const bool exitIfLeakPath,
                const refPtr<coordSetWriter>& leakPathFormatter,
                refPtr<surfaceWriter>& surfFormatter
            );

            //- Merge free-standing baffles
            void mergeFreeStandingBaffles
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
            );

            //- Split off (with optional buffer layers) unreachable areas
            //  of mesh. Does not introduce baffles.
            autoPtr<mapPolyMesh> splitMesh
            (
                const label nBufferLayers,
                const label nErodeCellZones,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,

                const pointField& locationsInMesh,
                const wordList& regionsInMesh,
                const pointField& locationsOutsideMesh,
                const bool exitIfLeakPath,
                const refPtr<coordSetWriter>& leakPathFormatter,
                refPtr<surfaceWriter>& surfFormatter
            );

            //- Remove cells from limitRegions if level -1
            autoPtr<mapPolyMesh> removeLimitShells
            (
                const label nBufferLayers,
                const label nErodeCellZones,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const pointField& locationsInMesh,
                const wordList& regionsInMesh,
                const pointField& locationsOutsideMesh
            );

            //- Find boundary points that connect to more than one cell
            //  region and split them.
            autoPtr<mapPolyMesh> dupNonManifoldPoints(const localPointRegion&);

            //- Find boundary points that connect to more than one cell
            //  region and split them.
            autoPtr<mapPolyMesh> dupNonManifoldPoints();

            //- Find boundary points that are on faceZones of type boundary
            //  and duplicate them
            autoPtr<mapPolyMesh> dupNonManifoldBoundaryPoints();

            //- Merge duplicate points
            autoPtr<mapPolyMesh> mergePoints
            (
                const labelList& pointToDuplicate
            );

            //- Create baffle for every internal face where ownPatch != -1.
            //  External faces get repatched according to ownPatch (neiPatch
            //  should be -1 for these)
            autoPtr<mapPolyMesh> createBaffles
            (
                const labelList& ownPatch,
                const labelList& neiPatch
            );

            //- Get zones of given type
            labelList getZones
            (
                const List<surfaceZonesInfo::faceZoneType>& fzTypes
            ) const;

            //- Subset baffles according to zones
            static List<labelPair> subsetBaffles
            (
                const polyMesh& mesh,
                const labelList& zoneIDs,
                const List<labelPair>& baffles
            );

            //- Map baffles after layer addition. Gets new-to-old face map.
            static void mapBaffles
            (
                const polyMesh& mesh,
                const labelList& faceMap,
                List<labelPair>& baffles
            );

            //- Get per-face information (faceZone, master/slave patch)
            void getZoneFaces
            (
                const labelList& zoneIDs,
                labelList& faceZoneID,
                labelList& ownPatch,
                labelList& neiPatch,
                labelList& nBaffles
            ) const;

            //- Create baffles for faces on faceZones. Return created baffles
            //  (= pairs of faces) and corresponding faceZone
            autoPtr<mapPolyMesh> createZoneBaffles
            (
                const labelList& zoneIDs,
                List<labelPair>& baffles,
                labelList& originatingFaceZone
            );

            //- Merge baffles. Gets pairs of faces and boundary faces to move
            //  onto (coupled) patches
            autoPtr<mapPolyMesh> mergeBaffles
            (
                const List<labelPair>&,
                const Map<label>& faceToPatch
            );

            //- Merge all baffles on faceZones
            autoPtr<mapPolyMesh> mergeZoneBaffles
            (
                const bool doInternalZones,
                const bool doBaffleZones
            );

            //- Put faces/cells into zones according to surface specification.
            //  Returns null if no zone surfaces present. Regions containing
            //  locationsInMesh/regionsInMesh will be put in corresponding
            //  cellZone. keepPoints is for backwards compatibility and sets
            //  all yet unassigned cells to be non-zoned (zone = -1)
            autoPtr<mapPolyMesh> zonify
            (
                const bool allowFreeStandingZoneFaces,
                const label nErodeCellZones,
                const pointField& locationsInMesh,
                const wordList& regionsInMesh,
                const pointField& locationsOutsideMesh,
                const bool exitIfLeakPath,
                const refPtr<coordSetWriter>& leakPathFormatter,
                refPtr<surfaceWriter>& surfFormatter,
                wordPairHashTable& zonesToFaceZone
            );


        // Other topo changes

            //- Helper:append patch to end of mesh.
            static label appendPatch
            (
                fvMesh&,
                const label insertPatchi,
                const word&,
                const dictionary&
            );

            //- Helper:add patch to mesh. Update all registered fields.
            //  Used by addMeshedPatch to add patches originating from surfaces.
            static label addPatch(fvMesh&, const word& name, const dictionary&);

            //- Add patch originating from meshing. Update meshedPatches_.
            label addMeshedPatch(const word& name, const dictionary&);

            //- Get patchIDs for patches added in addMeshedPatch.
            labelList meshedPatches() const;

            //- Add/lookup faceZone and update information. Return index of
            //  faceZone
            label addFaceZone
            (
                const word& fzName,
                const word& masterPatch,
                const word& slavePatch,
                const surfaceZonesInfo::faceZoneType& fzType
            );

            //- Lookup faceZone information. Return false if no information
            //  for faceZone
            bool getFaceZoneInfo
            (
                const word& fzName,
                label& masterPatchID,
                label& slavePatchID,
                surfaceZonesInfo::faceZoneType& fzType
            ) const;

            //- Add pointZone if does not exist. Return index of zone
            label addPointZone(const word& name);

            //- Count number of faces per patch edge. Parallel consistent.
            labelList countEdgeFaces(const uindirectPrimitivePatch& pp) const;

            //- Select coupled faces that are not collocated
            void selectSeparatedCoupledFaces(boolList&) const;

            //- Find any intersection of surface. Store in surfaceIndex_.
            void updateIntersections(const labelList& changedFaces);

            //- Calculate nearest intersection for selected mesh faces
            void nearestIntersection
            (
                const labelList& surfacesToTest,
                const labelList& testFaces,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2
            ) const;

            //- Remove cells. Put exposedFaces into exposedPatchIDs.
            autoPtr<mapPolyMesh> doRemoveCells
            (
                const labelList& cellsToRemove,
                const labelList& exposedFaces,
                const labelList& exposedPatchIDs,
                removeCells& cellRemover
            );

            //- Find cell point is in. Uses optional perturbation to re-test.
            //  Returns -1 on processors that do not have the cell.
            static label findCell
            (
                const polyMesh&,
                const vector& perturbVec,
                const point& p
            );

            //- Find region point is in. Uses optional perturbation to re-test.
            static label findRegion
            (
                const polyMesh&,
                const labelList& cellRegion,
                const vector& perturbVec,
                const point& p
            );

            //- Find regions points are in.
            //  \return number of cells to be removed
            static label findRegions
            (
                const polyMesh&,
                const vector& perturbVec,
                const pointField& locationsInMesh,
                const pointField& locationsOutsideMesh,
                const label nRegions,
                labelList& cellRegion,
                const boolList& blockedFace,
                // Leak-path
                const bool exitIfLeakPath,
                const refPtr<coordSetWriter>& leakPathFormatter
            );

            //- Split mesh. Keep part containing point. Return empty map if
            //  no cells removed.
            autoPtr<mapPolyMesh> splitMeshRegions
            (
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const pointField& locationsInMesh,
                const pointField& locationsOutsideMesh,
                // Leak-path
                const bool exitIfLeakPath,
                const refPtr<coordSetWriter>& leakPathFormatter
            );

            //- Helper: split face into:
            //          - f0 : split[0] to split[1]
            //          - f1 : split[1] to split[0]
            static void splitFace
            (
                const face& f,
                const labelPair& split,

                face& f0,
                face& f1
            );

            //- Split faces into two
            void doSplitFaces
            (
                const labelList& splitFaces,
                const labelPairList& splits,
                const labelPairList& splitPatches,
                polyTopoChange& meshMod
            ) const;

            //- Split faces along diagonal. Maintain mesh quality. Return
            //  total number of faces split.
            label splitFacesUndo
            (
                const labelList& splitFaces,
                const labelPairList& splits,
                const labelPairList& splitPatches,
                const dictionary& motionDict,

                labelList& duplicateFace,
                List<labelPair>& baffles
            );

            //- Update local numbering for mesh redistribution
            void distribute(const mapDistributePolyMesh&);

            //- Update for external change to mesh. changedFaces are in new mesh
            //  face labels.
            void updateMesh
            (
                const mapPolyMesh&,
                const labelList& changedFaces
            );

            //- Helper: reorder list according to map.
            template<class T>
            static void updateList
            (
                const labelList& newToOld,
                const T& nullValue,
                List<T>& elems
            );


            // Restoring : is where other processes delete and reinsert data.

                //- Signal points/face/cells for which to store data
                void storeData
                (
                    const labelList& pointsToStore,
                    const labelList& facesToStore,
                    const labelList& cellsToStore
                );

                //- Update local numbering + undo
                //  Data to restore given as new pointlabel + stored pointlabel
                //  (i.e. what was in pointsToStore)
                void updateMesh
                (
                    const mapPolyMesh&,
                    const labelList& changedFaces,
                    const Map<label>& pointsToRestore,
                    const Map<label>& facesToRestore,
                    const Map<label>& cellsToRestore
                );

            // Merging coplanar faces and edges

                //- Merge coplanar faces if sets are of size mergeSize
                //  (usually 4)
                label mergePatchFaces
                (
                    const scalar minCos,
                    const scalar concaveCos,
                    const label mergeSize,
                    const labelList& patchIDs,
                    const meshRefinement::FaceMergeType mergeType
                );

                //- Merge coplanar faces. preserveFaces is != -1 for faces
                //  to be preserved
                label mergePatchFacesUndo
                (
                    const scalar minCos,
                    const scalar concaveCos,
                    const labelList& patchIDs,
                    const dictionary& motionDict,
                    const labelList& preserveFaces,
                    const meshRefinement::FaceMergeType mergeType
                );

                autoPtr<mapPolyMesh> doRemovePoints
                (
                    removePoints& pointRemover,
                    const boolList& pointCanBeDeleted
                );

                autoPtr<mapPolyMesh> doRestorePoints
                (
                    removePoints& pointRemover,
                    const labelList& facesToRestore
                );

                labelList collectFaces
                (
                    const labelList& candidateFaces,
                    const labelHashSet& set
                ) const;

                // Pick up faces of cells of faces in set.
                labelList growFaceCellFace
                (
                    const labelUList& set
                ) const;

                // Pick up faces of cells of faces in set.
                labelList growFaceCellFace
                (
                    const labelHashSet& set
                ) const;

                //- Merge edges, maintain mesh quality. Return global number
                //  of edges merged
                label mergeEdgesUndo
                (
                    const scalar minCos,
                    const dictionary& motionDict
                );


        // Debug/IO

            //- Debugging: check that all faces still obey start()>end()
            void checkData();

            static void testSyncPointList
            (
                const string& msg,
                const polyMesh& mesh,
                const List<scalar>& fld
            );

            static void testSyncPointList
            (
                const string& msg,
                const polyMesh& mesh,
                const List<point>& fld
            );

            //- Compare two lists over all boundary faces
            template<class T>
            void testSyncBoundaryFaceList
            (
                const scalar mergeDistance,
                const string&,
                const UList<T>&,
                const UList<T>&
            ) const;

            //- Print list according to (collected and) sorted coordinate
            template<class T>
            static void collectAndPrint
            (
                const UList<point>& points,
                const UList<T>& data
            );

            //- Determine master point for subset of points. If coupled
            //  chooses only one
            static bitSet getMasterPoints
            (
                const polyMesh& mesh,
                const labelList& meshPoints
            );

            //- Determine master edge for subset of edges. If coupled
            //  chooses only one
            static bitSet getMasterEdges
            (
                const polyMesh& mesh,
                const labelList& meshEdges
            );

            //- Print some mesh stats.
            void printMeshInfo
            (
                const bool debug,
                const string& msg,
                const bool printCellLevel
            ) const;

            //- Replacement for Time::timeName() that returns oldInstance
            //- (if overwrite_)
            word timeName() const;

            //- Set instance of all local IOobjects
            void setInstance(const fileName&);

            //- Write mesh and all data
            bool write() const;

            //- Write refinement level as volScalarFields for postprocessing
            void dumpRefinementLevel() const;

            //- Debug: Write intersection information to OBJ format
            void dumpIntersections(const fileName& prefix) const;

            //- Do any one of above IO functions
            void write
            (
                const debugType debugFlags,
                const writeType writeFlags,
                const fileName&
            ) const;

            //- Helper: remove all relevant files from mesh instance
            static void removeFiles(const polyMesh&);

            //- Helper: calculate average
            template<class T>
            static T gAverage
            (
                const bitSet& isMasterElem,
                const UList<T>& values
            );

            //- Get/set write level
            static writeType writeLevel();
            static void writeLevel(const writeType);

            ////- Get/set output level
            //static outputType outputLevel();
            //static void outputLevel(const outputType);


            //- Helper: convert wordList into bit pattern using provided Enum
            template<class EnumContainer>
            static int readFlags
            (
                const EnumContainer& namedEnum,
                const wordList& words
            );

            //- Wrapper around dictionary::get which does not exit
            template<class Type>
            static Type get
            (
                const dictionary& dict,
                const word& keyword,
                const bool noExit,
                enum keyType::option matchOpt = keyType::REGEX,
                const Type& deflt = Zero
            );

            //- Wrapper around dictionary::subDict which does not exit
            static const dictionary& subDict
            (
                const dictionary& dict,
                const word& keyword,
                const bool noExit,
                enum keyType::option matchOpt = keyType::REGEX
            );

            //- Wrapper around dictionary::lookup which does not exit
            static ITstream& lookup
            (
                const dictionary& dict,
                const word& keyword,
                const bool noExit,
                enum keyType::option matchOpt = keyType::REGEX
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshRefinementTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
