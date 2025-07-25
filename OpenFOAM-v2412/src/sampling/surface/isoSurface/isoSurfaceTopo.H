/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::isoSurfaceTopo

Description
    Marching tet iso surface algorithm with optional filtering to keep only
    points originating from mesh edges.

SourceFiles
    isoSurfaceTopo.C
    isoSurfaceTopoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceTopo_H
#define isoSurfaceTopo_H

#include "bitSet.H"
#include "edgeList.H"
#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellShape;
class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceTopo Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceTopo
:
    public isoSurfaceBase
{
    // Private Data

        //- Per point: originating mesh vertex/cell centre combination.
        //  Vertices less than nPoints are mesh vertices,
        //  duplicate vertices indicate a snapped point
        edgeList pointToVerts_;


    // Private Classes

        //- Handling, bookkeeping for tet cuts
        class tetCutAddressing
        {
            // Bookkeeping hashes used during construction
            EdgeMap<label> vertsToPointLookup_;
            Map<label> snapVertsLookup_;

            // Filter information for face diagonals
            DynamicList<label> pointToFace_;
            DynamicList<bool> pointFromDiag_;

            // Final output
            DynamicList<edge> pointToVerts_;
            DynamicList<label> cutPoints_;

            //- List of cut (decomposed) cell tets. Debug use only.
            DynamicList<cellShape> debugCutTets_;

            bool debugCutTetsOn_;


        public:

        // Constructors

            //- Construct with reserved sizes
            tetCutAddressing
            (
                const label nCutCells,
                const bool useSnap,
                const bool useDebugCuts = false
            );


        // Member Functions

            //- Effective number of faces
            label nFaces() const { return cutPoints_.size()/3; }

            DynamicList<label>& pointToFace() { return pointToFace_; }
            DynamicList<bool>& pointFromDiag() { return pointFromDiag_; }

            DynamicList<edge>& pointToVerts() { return pointToVerts_; }
            DynamicList<label>& cutPoints() { return cutPoints_; }
            DynamicList<cellShape>& debugCutTets() { return debugCutTets_; }

            //- Number of debug cut tets
            label nDebugTets() const { return debugCutTets_.size(); }

            //- Debugging cut tets active
            bool debugCutTetsOn() const { return debugCutTetsOn_; }

            void clearDebug();
            void clearDiagonal();
            void clearHashes();

            //- Generate single point on edge
            label generatePoint
            (
                label facei,          //!< Originating mesh face for cut-point
                bool edgeIsDiagonal,  //!< Edge on face diagonal

                // 0: no snap, 1: snap first, 2: snap second
                const int snapEnd,

                const edge& vertices
            );

            //- Generate triangles from tet cut
            bool generatePoints
            (
                const label facei,
                const int tetCutIndex,  //!< Encoded tet cut + tet snapping
                const tetCell& tetLabels,
                const FixedList<bool, 6>& edgeIsDiagonal
            );
        };


    // Private Member Functions

        //- Generate triangle points from cell
        void generateTriPoints
        (
            const label celli,
            const bool isTet,
            const labelList& tetBasePtIs,

            tetCutAddressing& tetCutAddr
        ) const;


        // Simplification

            static void triangulateOutside
            (
                const bool filterDiag,
                const primitivePatch& pp,

                const boolUList& pointFromDiag,
                const labelUList& pointToFace,
                const label cellID,

                DynamicList<face>& compactFaces,
                DynamicList<label>& compactCellIDs
            );

            static void removeInsidePoints
            (
                Mesh& s,                            // Modify inplace
                const bool filterDiag,

                // Inputs
                const boolUList& pointFromDiag,
                const labelUList& pointToFace,
                const labelUList& start,            // Per cell:starting tri

                // Outputs
                DynamicList<label>& compactCellIDs  // Per face the cellID
            );


protected:

    // Editing

        //- Subset the surface using the selected faces.
        //
        //  \param[in] include the faces to select
        void inplaceSubsetMesh(const bitSet& include);


    // Sampling

        //- Interpolates cellData and pointData fields
        template<class Type>
        tmp<Field<Type>> interpolateTemplate
        (
            const Field<Type>& cellData,
            const Field<Type>& pointData
        ) const;

public:

    //- Runtime type information
    TypeName("isoSurfaceTopo");


    // Constructors

        //- Construct from cell and point values
        //
        //  \param ignoreCells cells to ignore in the cellValues
        //
        //  Control parameters include
        //  - bounds optional bounding box for trimming
        //  - mergeTol fraction of mesh bounding box for merging points
        isoSurfaceTopo
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const isoSurfaceParams& params = isoSurfaceParams(),
            const bitSet& ignoreCells = bitSet()
        );


    //- Destructor
    virtual ~isoSurfaceTopo() = default;


    // Member Functions

        //- Per point: originating mesh vertex/cell centre combination.
        //  Vertices less than nPoints are mesh vertices (encoding above),
        //  duplicate vertices indicate a snapped point
        const edgeList& pointToVerts() const noexcept
        {
            return pointToVerts_;
        }


    // Sampling

        declareIsoSurfaceInterpolateMethod(scalar);
        declareIsoSurfaceInterpolateMethod(vector);
        declareIsoSurfaceInterpolateMethod(sphericalTensor);
        declareIsoSurfaceInterpolateMethod(symmTensor);
        declareIsoSurfaceInterpolateMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isoSurfaceTopoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
