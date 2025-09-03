/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::snappyRefineDriver

Description

SourceFiles
    snappyRefineDriver.C

\*---------------------------------------------------------------------------*/

#ifndef snappyRefineDriver_H
#define snappyRefineDriver_H

#include "wordPairHashes.H"
#include "labelList.H"
#include "scalarField.H"
#include "Tuple2.H"
#include "coordSetWriter.H"
#include "DynamicList.H"
#include "labelVector.H"
#include "meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class refinementParameters;
class snapParameters;

class decompositionMethod;
class fvMeshDistribute;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class snappyRefineDriver Declaration
\*---------------------------------------------------------------------------*/

class snappyRefineDriver
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        fvMeshDistribute& distributor_;

        //- From surface region to patch
        const labelList globalToMasterPatch_;

        //- From surface region to patch
        const labelList globalToSlavePatch_;

        //- How to write lines. Used e.g. when writing leak-paths
        coordSetWriter& setFormatter_;

        //- How to write surfaces. Used e.g. when writing leak closure faces
        refPtr<surfaceWriter>& surfFormatter_;

        //- Are we operating in test mode?
        const bool dryRun_;


    // Private Member Functions

        void addNeighbours
        (
            const labelVector& n,
            const labelList& cellLevel,
            const labelVector& voxel,
            const label voxeli,
            DynamicList<labelVector>& front
        ) const;

        //- Rough estimate of cell size and cell count
        void estimateCellSizeAndCount
        (
            const refinementParameters& refineParams
        ) const;

        //- Refine all cells pierced by explicit feature edges
        label featureEdgeRefine
        (
            const refinementParameters& refineParams,
            const label maxIter,
            const label minRefine
        );

        //- Refine all cells containing small surface features
        label smallFeatureRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine all cells interacting with the surface
        label surfaceOnlyRefine
        (
            const refinementParameters& refineParams,
            const label maxIter,
            const label leakBlockageIter    // when to start leak closing
        );

        //- Refine all cells in small gaps
        label gapOnlyRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Block all cells in small gaps
        label surfaceProximityBlock
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine all cells in large gaps
        label bigGapOnlyRefine
        (
            const refinementParameters& refineParams,
            const bool spreadGapSize,
            const label maxIter
        );

        //- Refine cells with almost all sides refined
        label danglingCellRefine
        (
            const refinementParameters& refineParams,
            const label nFaces,
            const label maxIter
        );

        //- Refine cells with opposite faces with differing refinement level
        label refinementInterfaceRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Helper: see if any element in f has higher level than cLevel
        bool usesHigherLevel
        (
            const labelUList& boundaryPointLevel,
            const labelUList& f,
            const label cLevel
        ) const;

        //- Refine cells with a point/edge but not face on the boundary
        label boundaryRefinementInterfaceRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Remove all cells within intersected region
        void removeInsideCells
        (
            const refinementParameters& refineParams,
            const label nBufferLayers
        );

        //- Refine all cells inside/outside shell
        label shellRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        // Directional refinement and smoothing

            //- Refine (directional) all cells inside/outside shell
            label directionalShellRefine
            (
                const refinementParameters& refineParams,
                const label maxIter
            );

            //- Calculate local edge length from cell volumes
            void mergeAndSmoothRatio
            (
                const scalarList& allSeedPointDist,
                const label nSmoothExpansion,
                List<Tuple2<scalar, scalar>>&  keyAndValue
            );

            //- Smooth the directional expansion ratio
            label directionalSmooth(const refinementParameters& refineParams);


        //- Add baffles and remove unreachable cells
        void baffleAndSplitMesh
        (
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Add zones
        void zonify
        (
            const refinementParameters& refineParams,
            wordPairHashTable& zonesToFaceZone
        );

        void splitAndMergeBaffles
        (
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Erode sticking-out zone faces. Bit heuristic.
        void erodeNonManifoldZoneFaces
        (
            const refinementParameters& refineParams
        );

        //- Merge refined boundary faces (from exposing coarser cell)
        void mergePatchFaces
        (
            const meshRefinement::FaceMergeType mergeType,
            const refinementParameters& refineParams,
            const dictionary& motionDict
        );

        //- Optionally delete some small regions
        void deleteSmallRegions(const refinementParameters&);


        //- No copy construct
        snappyRefineDriver(const snappyRefineDriver&) = delete;

        //- No copy assignment
        void operator=(const snappyRefineDriver&) = delete;


public:

    //- Runtime type information
    ClassName("snappyRefineDriver");


    // Constructors

        //- Construct from components
        snappyRefineDriver
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor,
            const labelUList& globalToMasterPatch,
            const labelUList& globalToSlavePatch,
            coordSetWriter& setFormatter,
            refPtr<surfaceWriter>& surfFormatter,
            const bool dryRun = false
        );


    // Member Functions

        //- Do all the refinement
        void doRefine
        (
            const dictionary& refineDict,
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool prepareForSnapping,
            const meshRefinement::FaceMergeType mergeType,
            const dictionary& motionDict
        );

        //- Helper: add faceZones and patches
        static void addFaceZones
        (
            meshRefinement& meshRefiner,
            const refinementParameters& refineParams,
            const HashTable<Pair<word>>& faceZoneToPatches
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
