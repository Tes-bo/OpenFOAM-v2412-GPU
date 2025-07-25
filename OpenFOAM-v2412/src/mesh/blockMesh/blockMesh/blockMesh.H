/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::blockMesh

Description
    A multi-block mesh generator

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        prescale    | Point scaling before transform        | no  | 1.0
        scale       | Point scaling after transform         | no  | 1.0
        transform   | Point transform (origin, rotation)    | no  |
        vertices    |                                       | yes |
        blocks      |                                       | yes |
        edges       |                                       | no  |
        faces       |                                       | no  |
        boundary    | Boundary definition                   | no  |
        patches     | Alternate version for "boundary"      | no  |
        namedBlocks |                                       | no  |
        namedVertices |                                     | no  |
        mergeType   | Merging "points" or "topology"        | no  | topology
        checkFaceCorrespondence |                           | no  | true
        verbose     |                                       | no  | true
    \endtable

Note
    The \c prescale and \c scale can be a single scalar or a vector of
    values.

    The vertices, cells and patches for filling the blocks are demand-driven.

SourceFiles
    blockMesh.C
    blockMeshCheck.C
    blockMeshCreate.C
    blockMeshMerge.C
    blockMeshTopology.C

\*---------------------------------------------------------------------------*/

#ifndef blockMesh_H
#define blockMesh_H

#include "Enum.H"
#include "Switch.H"
#include "block.H"
#include "PtrList.H"
#include "cartesianCS.H"
#include "searchableSurfaces.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "blockVertexList.H"
#include "blockEdgeList.H"
#include "blockFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blockMesh Declaration
\*---------------------------------------------------------------------------*/

class blockMesh
:
    public PtrList<block>
{
public:

    // Typedefs

        //- The list of blocks is stored as a PtrList
        typedef PtrList<block> blockList;


    // Data Types

        //- The block merging strategy
        enum mergeStrategy
        {
            DEFAULT_MERGE,  //!< Default (TOPOLOGY), not selectable
            MERGE_TOPOLOGY, //!< "topology" merge by block topology (default)
            MERGE_POINTS    //!< "points" merge by point geometry
        };


private:

    // Static Data

        //- Names corresponding to the mergeStrategy
        static const Enum<mergeStrategy> strategyNames_;


    // Data Types

        //- Point transformations. Internal usage
        enum transformTypes : unsigned
        {
            NO_TRANSFORM = 0,   //!< No transformations
            ROTATION = 0x1,     //!< Local coordinate system rotation
            TRANSLATION = 0x2,  //!< Local coordinate system translation
            PRESCALING = 0x4,   //!< Uniform scale before rotate/translate
            PRESCALING3 = 0x8,  //!< Non-uniform scale before rotate/translate
            SCALING = 0x10,     //!< Uniform scale after rotate/translate
            SCALING3 = 0x20     //!< Non-uniform scale after rotate/translate
        };


    // Private Data

        //- Reference to mesh dictionary
        const IOdictionary& meshDict_;

        //- Output verbosity level
        int verbose_;

        //- Check face consistency (defaults to true)
        bool checkFaceCorrespondence_;

        //- The mesh merging strategy
        mergeStrategy mergeStrategy_;

        //- Types of point transformations requested
        unsigned transformType_;

        //- Optional searchable geometry to project face-points to
        searchableSurfaces geometry_;

        //- The list of block vertices
        blockVertexList blockVertices_;

        //- The list of block vertex positions
        pointField vertices_;

        //- The list of curved edges
        blockEdgeList edges_;

        //- The list of curved faces
        blockFaceList faces_;

        //- The scaling factor before rotate/translate
        vector prescaling_;

        //- The scaling factor after rotate/translate
        vector scaling_;

        //- Local coordinate system transformation
        coordSystem::cartesian transform_;

        //- The blocks themselves (the topology) as a polyMesh
        autoPtr<polyMesh> topologyPtr_;

        //- The sum of all cells in each block
        label nPoints_;

        //- The sum of all cells in each block
        label nCells_;

        //- The merge points information
        labelList mergeList_;

        //- The point offset added to each block. Offsets into mergeList_
        labelList blockOffsets_;

        mutable pointField points_;

        mutable cellShapeList cells_;

        mutable faceListList patches_;


    // Private Member Functions

        //- Get scaling and/or coordinate transforms
        //  \return True if scaling and/or transformations are needed
        bool readPointTransforms(const dictionary& dict);

        void readPatches
        (
            const dictionary& meshDescription,
            faceListList& tmpBlocksPatches,
            wordList& patchNames,
            wordList& patchTypes,
            wordList& nbrPatchNames
        );

        void readBoundary
        (
            const dictionary& meshDescription,
            wordList& patchNames,
            faceListList& tmpBlocksPatches,
            PtrList<dictionary>& patchDicts
        );

        //- Topology blocks as cell shapes
        cellShapeList getBlockShapes() const;

        autoPtr<polyMesh> createTopology
        (
            const IOdictionary& dict,
            const word& regionName
        );


        //- Simple checks for collapsed hex cells
        bool checkDegenerate() const;

        void check(const polyMesh& bm, const dictionary& dict) const;

        //- Determine merge info and final number of cells/points
        //- based on point distances
        void calcGeometricalMerge();

        //- Determine merge info and final number of cells/points
        //- based on block topology
        void calcTopologicalMerge();

        faceList createPatchFaces(const polyPatch& patchTopologyFaces) const;

        void createPoints() const;
        void createCells() const;
        void createPatches() const;


        //- No copy construct
        blockMesh(const blockMesh&) = delete;

        //- No copy assignment
        void operator=(const blockMesh&) = delete;


public:

    // Static Data

        //- The default verbosity (true)
        static bool verboseOutput;


    //- Runtime type information
    ClassName("blockMesh");


    // Constructors

        //- Construct from IOdictionary for given region
        //  Default is topological merging.
        explicit blockMesh
        (
            const IOdictionary& dict,
            const word& regionName = polyMesh::defaultRegion,
            mergeStrategy strategy = mergeStrategy::DEFAULT_MERGE,
            int verbosity = 0   // 0: use static or dictionary value
        );


    //- Destructor
    ~blockMesh() = default;


    // Member Functions

    // General Access, Description

        //- Access to input dictionary
        const dictionary& meshDict() const noexcept
        {
            return meshDict_;
        }

        //- Optional searchable geometry to project face-points to
        const searchableSurfaces& geometry() const noexcept
        {
            return geometry_;
        }

        //- The curved edges
        const blockEdgeList& edges() const noexcept
        {
            return edges_;
        }

        //- The curved faces
        const blockFaceList& faces() const noexcept
        {
            return faces_;
        }

        //- True if the blockMesh topology exists
        bool good() const noexcept;

        //- Same as good()
        bool valid() const noexcept { return good(); }

        //- Return the patch names
        wordList patchNames() const;

        //- Number of blocks with specified zones
        label numZonedBlocks() const;


    // Point Transformations

        //- True if scaling and/or transformations are needed
        bool hasPointTransforms() const noexcept;

        //- Apply coordinate transforms and scaling
        bool inplacePointTransforms(pointField& pts) const;

        //- Apply coordinate transforms and scaling
        tmp<pointField> globalPosition(const pointField& localPoints) const;


    // Block Topology

        //- Reference to point field defining the blocks,
        //- these points are \b unscaled and \b non-transformed
        const pointField& vertices() const noexcept;

        //- Point field defining the blocks,
        //- optionally transformed and scaled
        tmp<pointField> vertices(bool applyTransform) const;

        //- The blockMesh topology as a polyMesh
        //- \b unscaled and \b non-transformed
        const polyMesh& topology() const;

        //- The blockMesh topology as a polyMesh
        //- optionally transformed and scaled
        refPtr<polyMesh> topology(bool applyTransform) const;


    // Detailed Mesh

        //- The points for the entire mesh.
        //- These points \b are scaled and transformed
        const pointField& points() const;

        //- Return cell shapes list
        const cellShapeList& cells() const;

        //- Return the patch face lists
        const faceListList& patches() const;

        //- Patch information from the topology mesh
        PtrList<dictionary> patchDicts() const;


    // Verbosity

        //- Output verbosity level
        int verbose() const noexcept;

        //- Change the output verbosity level.
        //  \return old level
        int verbose(const int level) noexcept;


    // Mesh Generation

        //- Create polyMesh, with cell zones
        autoPtr<polyMesh> mesh(const IOobject& io) const;


    // Housekeeping

        //- Old (v2106 and earlier) uniform scaling factor
        //  \deprecated use inplacePointTransforms or globalPosition instead
        scalar scaleFactor() const
        {
            return scaling_.x();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
