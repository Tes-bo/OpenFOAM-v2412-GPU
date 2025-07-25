/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::snappySnapDriver

Description
    All to do with snapping to surface

SourceFiles
    snappySnapDriver.C
    snappySnapDriverFeature.C

\*---------------------------------------------------------------------------*/

#ifndef snappySnapDriver_H
#define snappySnapDriver_H

#include "meshRefinement.H"
#include "DynamicField.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class motionSmoother;
class refinementParameters;
class snapParameters;
class pointConstraint;
class layerParameters;
class displacementMotionSolver;
class addPatchCellLayer;

/*---------------------------------------------------------------------------*\
                           Class snappySnapDriver Declaration
\*---------------------------------------------------------------------------*/

class snappySnapDriver
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- From global surface region to master side patch
        const labelList globalToMasterPatch_;

        //- From global surface region to slave side patch
        const labelList globalToSlavePatch_;

        //- Are we operating in test mode?
        const bool dryRun_;


    // Private Member Functions


        // Snapping

            //- Top-level: snap onto surface & features and add buffer layer
            void doSnapBufferLayers
            (
                const dictionary& snapDict,
                const dictionary& motionDict,
                const meshRefinement::FaceMergeType mergeType,
                const scalar featureCos,
                const scalar planarAngle,
                const snapParameters& snapParams
            );

            //- Calculates (geometric) shared points
            //  Requires bitSet to be sized and initialised
            static label getCollocatedPoints
            (
                const scalar tol,
                const pointField&,
                bitSet&
            );

            //- Calculate displacement (over all mesh points) to move points
            //  to average of connected cell centres
            static tmp<pointField> smoothInternalDisplacement
            (
                const meshRefinement& meshRefiner,
                const motionSmoother&
            );

            //- Calculate displacement per patch point to smooth out patch.
            //  Quite complicated in determining which points to move where.
            static tmp<pointField> smoothPatchDisplacement
            (
                const motionSmoother&,
                const List<labelPair>&
            );

            static tmp<pointField> avg
            (
                const polyMesh& mesh,
                const bitSet& isMasterPoint,
                const indirectPrimitivePatch&,
                const pointField&
            );

            //- Calculate displacement per patch point. Wip.
            static tmp<pointField> smoothLambdaMuPatchDisplacement
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& pp,
                const List<labelPair>& baffles
            );

            tmp<scalarField> wantedThickness
            (
                const indirectPrimitivePatch& pp,
                const scalar cellSizeFraction
            ) const;

            //- Create/update pointMesh
            const pointMesh& makePointMesh
            (
                const indirectPrimitivePatch& pp,
                const pointConstraintList& pointConstraints,
                const word& allEdgePatchName,
                const word& allPointPatchName
            ) const;

            autoPtr<displacementMotionSolver> makeMotionSolver
            (
                const pointMesh& pMesh,
                const dictionary& snapDict,
                const labelList& adaptPatchIDs
            ) const;

            void setDisplacement
            (
                const indirectPrimitivePatch& pp,
                const pointField& patchDisp,
                const labelList& adaptPatchIDs,
                const pointField& points0,
                pointVectorField& fld
            );

            autoPtr<mapPolyMesh> addBufferLayers
            (
                const indirectPrimitivePatch& pp,

                //const scalar cellSizeFraction,
                const pointField& thickness,
                // Layer mesh modifier
                addPatchCellLayer& addLayer
            );

            //- Check that face zones are synced
            void checkCoupledFaceZones() const;

            //- Per edge distance to patch
            static tmp<scalarField> edgePatchDist
            (
                const pointMesh&,
                const indirectPrimitivePatch&
            );

            //- Write displacement as .obj file.
            static void dumpMove
            (
                const fileName&,
                const pointField&,
                const pointField&
            );

            //- Check displacement is outwards pointing
            static bool outwardsDisplacement
            (
                const indirectPrimitivePatch&,
                const vectorField&
            );

            //- Freeze points on pointZone or (inside of) faceZone
            static void freezeExposedPoints
            (
                const meshRefinement& meshRefiner,
                const word& fzName,     // faceZone name
                const word& pzName,     // pointZone name
                const indirectPrimitivePatch& outside,
                vectorField& patchDisp
            );

            //- Detect warpage
            void detectWarpedFaces
            (
                const scalar featureCos,
                const indirectPrimitivePatch& pp,

                DynamicList<label>& splitFaces,
                DynamicList<labelPair>& splits
            ) const;

            //- Get per face -1 or label of opposite face if on internal/baffle
            //  faceZone
            labelList getInternalOrBaffleDuplicateFace() const;

            //- Get points both on patch and facezone.
            static void getZoneSurfacePoints
            (
                const fvMesh& mesh,
                const indirectPrimitivePatch&,
                const word& zoneName,

                bitSet& pointOnZone
            );

            //- Get points both on patch and facezone.
            template<class FaceList>
            static labelList getFacePoints
            (
                const indirectPrimitivePatch& pp,
                const FaceList& faces
            );

            //- Per patch point calculate point on nearest surface.
            //  Return displacement of patch points.
            static void calcNearestSurface
            (
                const refinementSurfaces& surfaces,

                const labelList& surfacesToTest,
                const labelListList& regionsToTest,

                const pointField& localPoints,
                const labelList& zonePointIndices,

                scalarField& minSnapDist,
                labelList& snapSurf,
                vectorField& patchDisp,

                // Optional: nearest point, normal
                pointField& nearestPoint,
                vectorField& nearestNormal
            );


            // Feature line snapping

                //- Is point on two feature edges that make a largish angle?
                bool isFeaturePoint
                (
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const bitSet& isFeatureEdge,
                    const label pointi
                ) const;

                void smoothAndConstrain
                (
                    const bitSet& isMasterEdge,
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const List<pointConstraint>& constraints,
                    vectorField& disp
                ) const;
                //void smoothAndConstrain2
                //(
                //    const bool applyConstraints,
                //    const indirectPrimitivePatch& pp,
                //    const List<pointConstraint>& constraints,
                //    vectorField& disp
                //) const;
                void calcNearest
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    vectorField& pointDisp,
                    vectorField& pointSurfaceNormal,
                    vectorField& pointRotation
                ) const;
                void calcNearestFace
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& faceSnapDist,
                    vectorField& faceDisp,
                    vectorField& faceSurfaceNormal,
                    labelList& faceSurfaceRegion
                    //vectorField& faceRotation
                ) const;
                void calcNearestFacePointProperties
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,

                    const vectorField& faceDisp,
                    const vectorField& faceSurfaceNormal,
                    const labelList& faceSurfaceRegion,

                    List<List<point>>& pointFaceSurfNormals,
                    List<List<point>>& pointFaceDisp,
                    List<List<point>>& pointFaceCentres,
                    List<labelList>&    pointFacePatchID
                ) const;
                void correctAttraction
                (
                    const DynamicList<point>& surfacePoints,
                    const DynamicList<label>& surfaceCounts,
                    const point& edgePt,
                    const vector& edgeNormal,   // normalised normal
                    const point& pt,
                    vector& edgeOffset  // offset from pt to point on edge
                ) const;


                //- For any reverse (so from feature back to mesh) attraction:
                //  add attraction if diagonal points on face attracted
                void stringFeatureEdges
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Remove constraints of points next to multi-patch points
                //  to give a bit more freedom of the mesh to conform to the
                //  multi-patch points. Bit dodgy for simple cases.
                void releasePointsNextToMultiPatch
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,

                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Detect any diagonal attraction. Returns indices in face
                //  or (-1, -1) if none
                labelPair findDiagonalAttraction
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& patchAttraction,
                    const List<pointConstraint>& patchConstraints,
                    const label facei
                ) const;

                scalar pyrVol
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const face& localF,
                    const point& cc
                ) const;
                void facePoints
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const vectorField& nearestAttraction,
                    const face& f,
                    DynamicField<point>& points
                ) const;
                scalar pyrVol
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const vectorField& nearestAttraction,
                    const face& localF,
                    const point& cc
                ) const;
                Tuple2<point, vector> centreAndNormal
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const vectorField& nearestAttraction,
                    const face& localF
                ) const;
                bool isSplitAlignedWithFeature
                (
                    const scalar featureCos,
                    const point& newPt0,
                    const pointConstraint& pc0,
                    const point& newPt1,
                    const pointConstraint& pc1
                ) const;
                bool isConcave
                (
                    const point& c0,
                    const vector& area0,
                    const point& c1,
                    const vector& area1,
                    const scalar concaveCos
                ) const;
                labelPair findDiagonalAttraction
                (
                    const scalar featureCos,
                    const scalar concaveCos,
                    const scalar minAreaFraction,
                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const vectorField& patchAttraction,
                    const List<pointConstraint>& patchConstraints,
                    const vectorField& nearestAttraction,
                    const vectorField& nearestNormal,
                    const label faceI,

                    DynamicField<point>& points0,
                    DynamicField<point>& points1
                ) const;

                //- Do all logic on whether to add face cut to diagonal
                //  attraction
                void splitDiagonals
                (
                    const scalar featureCos,
                    const scalar concaveCos,
                    const scalar minAreaFraction,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const vectorField& nearestAttraction,
                    const vectorField& nearestNormal,
                    const List<labelList>& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints,
                    DynamicList<label>& splitFaces,
                    DynamicList<labelPair>& splits,
                    DynamicList<labelPair>& splitPatches
                ) const;

                //- Avoid attraction across face diagonal since would
                //  cause face squeeze
                void avoidDiagonalAttraction
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Write some stats about constraints
                void writeStats
                (
                    const indirectPrimitivePatch& pp,
                    const bitSet& isMasterPoint,
                    const List<pointConstraint>& patchConstraints
                ) const;

                //- Return hit if on multiple points
                pointIndexHit findMultiPatchPoint
                (
                    const point& pt,
                    const labelList& patchIDs,
                    const List<point>& faceCentres
                ) const;

                //- Return hit if faces-on-the-same-normalplane are on multiple
                //  patches
                pointIndexHit findMultiPatchPoint
                (
                    const point& pt,
                    const labelList& pfPatchID,
                    const DynamicList<vector>& surfaceNormals,
                    const labelList& faceToNormalBin
                ) const;

                //- Return index of similar normal
                label findNormal
                (
                    const scalar featureCos,
                    const vector& faceSurfaceNormal,
                    const DynamicList<vector>& surfaceNormals
                ) const;

                //- Determine attraction and constraints for single point
                //  using sampled surrounding of the point
                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    const label pointi,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    DynamicList<point>& surfacePoints,
                    DynamicList<vector>& surfaceNormals,
                    labelList& faceToNormalBin,

                    vector& patchAttraction,
                    pointConstraint& patchConstraint
                ) const;

                //- Determine attraction and constraints for all points
                //  using sampled surrounding of the point
                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Determine geometric features and attraction to equivalent
                //  surface features
                void determineFeatures
                (
                    const label iter,
                    const scalar featureCos,
                    const bool multiRegionFeatureSnap,
                    const bool strictRegionFeatureSnap, // special feat-point

                    const indirectPrimitivePatch&,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Determine features originating from bafles and
                //  and add attraction to equivalent surface features
                void determineBaffleFeatures
                (
                    const label iter,
                    const bool baffleFeaturePoints,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;
                void reverseAttractMeshPoints
                (
                    const label iter,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,

                    // Feature-point to pp point
                    const List<labelList>& pointAttractor,
                    const List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    const List<List<DynamicList<point>>>& edgeAttractors,
                    const List<List<DynamicList<pointConstraint>>>&,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Find point on nearest feature edge (within searchDist).
                //  Return point and feature
                //  and store feature-edge to mesh-point and vice versa
                Tuple2<label, pointIndexHit> findNearFeatureEdge
                (
                    const bool isRegionEdge,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,
                    const label pointi,
                    const point& estimatedPt,

                    List<List<DynamicList<point>>>&,
                    List<List<DynamicList<pointConstraint>>>&,
                    vectorField&,
                    List<pointConstraint>&
                ) const;

                //- Find nearest feature point (within searchDist).
                //  Return feature point
                //  and store feature-point to mesh-point and vice versa.
                //  If another mesh point already referring to this feature
                //  point and further away, reset that one to a near feature
                //  edge (using findNearFeatureEdge above)
                Tuple2<label, pointIndexHit> findNearFeaturePoint
                (
                    const bool isRegionEdge,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,
                    const label pointi,
                    const point& estimatedPt,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void featureAttractionUsingFeatureEdges
                (
                    const label iter,
                    const bool multiRegionFeatureSnap,
                    const bool strictRegionFeatureSnap,

                    const bool detectBaffles,
                    const bool baffleFeaturePoints,
                    const bool releasePoints,
                    const bool stringFeatures,
                    const bool avoidDiagonal,

                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    const vectorField& nearestNormal,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void preventFaceSqueeze
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    const scalarField& snapDist,
                    const vectorField& nearestAttraction,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Top level feature attraction routine. Gets given
                //  displacement to nearest surface in nearestDisp
                //  and calculates new displacement taking into account
                //  features
                vectorField calcNearestSurfaceFeature
                (
                    const snapParameters& snapParams,
                    const bool alignMeshEdges,
                    const bool strictRegionFeatureSnap, //points on >=3 patches
                    const label iter,
                    const scalar featureCos,
                    const scalar featureAttract,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    const vectorField& nearestNormal,
                    const indirectPrimitivePatch& pp,
                    const pointField& ppLocalPoints,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints,

                    DynamicList<label>& splitFaces,
                    DynamicList<labelPair>& splits,
                    DynamicList<labelPair>& splitPatches
                ) const;


        //- No copy construct
        snappySnapDriver(const snappySnapDriver&) = delete;

        //- No copy assignment
        void operator=(const snappySnapDriver&) = delete;


public:

    //- Runtime type information
    ClassName("snappySnapDriver");


    // Constructors

        //- Construct from components
        snappySnapDriver
        (
            meshRefinement& meshRefiner,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch,
            const bool dryRun = false
        );


    // Member Functions

        // Snapping

            //- Merge baffles.
            autoPtr<mapPolyMesh> mergeZoneBaffles(const List<labelPair>&);

            //- Calculate edge length per patch point.
            static scalarField calcSnapDistance
            (
                const fvMesh& mesh,
                const snapParameters& snapParams,
                const indirectPrimitivePatch&
            );

            //- Smooth the mesh (patch and internal) to increase visibility
            //  of surface points (on castellated mesh) w.r.t. surface.
            static void preSmoothPatch
            (
                const meshRefinement& meshRefiner,
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            //- Helper: calculate average cell centre per point
            static tmp<pointField> avgCellCentres
            (
                const fvMesh& mesh,
                const indirectPrimitivePatch&
            );

            //- Per patch point override displacement if in gap situation
            void detectNearSurfaces
            (
                const scalar planarCos,
                const indirectPrimitivePatch& pp,
                const pointField& ppLocalPoints,
                const pointField& nearestPoint,
                const vectorField& nearestNormal,
                vectorField& disp
            ) const;

            //- Per patch point calculate point on nearest surface. Set as
            //  boundary conditions of motionSmoother displacement field. Return
            //  displacement of patch points.
            static vectorField calcNearestSurface
            (
                const bool strictRegionSnap,
                const meshRefinement& meshRefiner,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const indirectPrimitivePatch& pp,
                const pointField& ppLocalPoints,
                const scalarField& snapDist,
                pointField& nearestPoint,
                vectorField& nearestNormal
            );

            ////- Per patch point calculate point on nearest surface. Set as
            ////  boundary conditions of motionSmoother displacement field.
            ////  Return displacement of patch points.
            //static vectorField calcNearestLocalSurface
            //(
            //    const meshRefinement& meshRefiner,
            //    const scalarField& snapDist,
            //    const indirectPrimitivePatch&
            //);

            //- Smooth the displacement field to the internal.
            void smoothDisplacement
            (
                const snapParameters& snapParams,
                motionSmoother&
            ) const;

            //- Do the hard work: move the mesh according to displacement,
            //  locally relax the displacement. Return true if ended up with
            //  correct mesh, false if not.
            bool scaleMesh
            (
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            //- Repatch faces according to surface nearest the face centre
            //  - calculate face-wise snap distance as max of point-wise
            //  - calculate face-wise nearest surface point
            //  - repatch face according to patch for surface point.
            autoPtr<mapPolyMesh> repatchToSurface
            (
                const snapParameters& snapParams,
                const labelList& adaptPatchIDs,
                const labelList& preserveFaces
            );

            //- Snap onto surface & features
            void doSnap
            (
                const dictionary& snapDict,
                const dictionary& motionDict,
                const meshRefinement::FaceMergeType mergeType,
                const scalar featureCos,
                const scalar planarAngle,
                const snapParameters& snapParams
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "snappySnapDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
