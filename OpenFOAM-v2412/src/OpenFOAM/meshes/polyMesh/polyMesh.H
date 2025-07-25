/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::polyMesh

Description
    Mesh consisting of general polyhedral cells.

SourceFiles
    polyMesh.C
    polyMeshInitMesh.C
    polyMeshClear.C
    polyMeshFromShapeMesh.C
    polyMeshIO.C
    polyMeshUpdate.C
    polyMeshCheck.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_polyMesh_H
#define Foam_polyMesh_H

#include "objectRegistry.H"
#include "primitiveMesh.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"
#include "cellShapeList.H"
#include "pointIOField.H"
#include "faceIOList.H"
#include "labelIOList.H"
#include "polyBoundaryMesh.H"
#include "boundBox.H"
#include "globalIndex.H"
#include "pointZoneMesh.H"
#include "faceZoneMesh.H"
#include "cellZoneMesh.H"
#include "meshState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class globalMeshData;
class mapPolyMesh;
class polyMeshTetDecomposition;
class treeDataCell;
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                          Class polyMesh Declaration
\*---------------------------------------------------------------------------*/

class polyMesh
:
    public objectRegistry,
    public primitiveMesh
{
public:

    // Public Data

        //- Enumeration defining the state of the mesh after a read update.
        //  Used for post-processing applications, where the mesh
        //  needs to update based on the files written in time
        //  directories
        enum readUpdateState
        {
            UNCHANGED,
            POINTS_MOVED,
            TOPO_CHANGE,
            TOPO_PATCH_CHANGE
        };

        //- Enumeration defining the decomposition of the cell for
        //  inside/outside test
        enum cellDecomposition
        {
            FACE_PLANES,      //- Faces considered as planes

            FACE_CENTRE_TRIS, //- Faces decomposed into triangles
                              //  using face-centre

            FACE_DIAG_TRIS,   //- Faces decomposed into triangles diagonally

            CELL_TETS         //- Cell decomposed into tets
        };


private:

    // Private Data

        //- Mesh and solver state data
        meshState data_;


        // Primitive mesh data

            //- Points
            pointIOField points_;

            //- Faces
            faceCompactIOList faces_;

            //- Face owner
            labelIOList owner_;

            //- Face neighbour
            labelIOList neighbour_;

            //- Have the primitives been cleared
            bool clearedPrimitives_;


            //- Boundary mesh
            mutable polyBoundaryMesh boundary_;

            //- Mesh bounding-box.
            //  Created from points on construction, updated when the mesh moves
            boundBox bounds_;

            //- Communicator used for parallel communication
            label comm_;

            //- Vector of non-constrained directions in mesh
            //  defined according to the presence of empty and wedge patches
            mutable Vector<label> geometricD_;

            //- Vector of valid directions in mesh
            //  defined according to the presence of empty patches
            mutable Vector<label> solutionD_;

            //- Base point for face decomposition into tets
            mutable autoPtr<labelIOList> tetBasePtIsPtr_;

            //- Search tree to allow spatial cell searching
            mutable autoPtr<indexedOctree<treeDataCell>> cellTreePtr_;


        // Zoning information

            //- Point zones
            pointZoneMesh pointZones_;

            //- Face zones
            faceZoneMesh faceZones_;

            //- Cell zones
            cellZoneMesh cellZones_;


        //- Parallel info
        mutable autoPtr<globalMeshData> globalMeshDataPtr_;


        // Mesh motion related data

            //- Is the mesh moving
            bool moving_;

            //- Is the mesh topology changing
            bool topoChanging_;

            //- Store old cell centres?
            mutable bool storeOldCellCentres_;

            //- Current time index for mesh motion
            mutable label curMotionTimeIndex_;

            //- Old points (for the last mesh motion)
            mutable autoPtr<pointField> oldPointsPtr_;

            //- Old cell centres (for the last mesh motion)
            mutable autoPtr<pointField> oldCellCentresPtr_;



    // Private Member Functions

        //- No copy construct
        polyMesh(const polyMesh&) = delete;

        //- No copy assignment
        void operator=(const polyMesh&) = delete;

        //- Initialise the polyMesh from the primitive data
        void initMesh();

        //- Initialise the polyMesh from the given set of cells
        void initMesh(cellList& c);

        //- Calculate the valid directions in the mesh from the boundaries
        void calcDirections() const;

        //- Calculate the cell shapes from the primitive
        //  polyhedral information
        void calcCellShapes() const;

        //- Read and return the tetBasePtIs
        autoPtr<labelIOList> readTetBasePtIs() const;


        // Helper functions for constructor from cell shapes

            labelListList cellShapePointCells(const cellShapeList&) const;

            labelList facePatchFaceCells
            (
                const faceList& patchFaces,
                const labelListList& pointCells,
                const faceListList& cellsFaceShapes,
                const label patchID
            ) const;

            void setTopology
            (
                const cellShapeList& cellsAsShapes,
                const faceListList& boundaryFaces,
                const wordList& boundaryPatchNames,
                labelList& patchSizes,
                labelList& patchStarts,
                label& defaultPatchStart,
                label& nFaces,
                cellList& cells
            );


        // Geometry checks

            //- Check non-orthogonality
            bool checkFaceOrthogonality
            (
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                const bool detailedReport,
                labelHashSet* setPtr
            ) const;

            //- Check face skewness
            bool checkFaceSkewness
            (
                const pointField& points,
                const vectorField& fCtrs,
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                const bool detailedReport,
                labelHashSet* setPtr
            ) const;

            bool checkEdgeAlignment
            (
                const pointField& p,
                const bool report,
                const Vector<label>& directions,
                labelHashSet* setPtr
            ) const;

            bool checkCellDeterminant
            (
                const vectorField& faceAreas,
                const bool report,
                labelHashSet* setPtr,
                const Vector<label>& meshD
            ) const;

            bool checkFaceWeight
            (
                const vectorField& fCtrs,
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                const scalar minWeight,
                labelHashSet* setPtr
            ) const;

            bool checkVolRatio
            (
                const scalarField& cellVols,
                const bool report,
                const scalar minRatio,
                labelHashSet* setPtr
            ) const;


public:

    // Public Typedefs

        //- The mesh type
        typedef polyMesh Mesh;

        //- The boundary type associated with the mesh
        typedef polyBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("polyMesh");

    //- Return the default region name
    static word defaultRegion;

    //- Return the mesh sub-directory name (usually "polyMesh")
    static word meshSubDir;


    // Constructors

        //- Read construct from IOobject
        explicit polyMesh(const IOobject& io, const bool doInit = true);

        //- Construct from IOobject or as zero-sized mesh
        //  Boundary is added using addPatches() member function
        polyMesh(const IOobject& io, const Foam::zero, const bool syncPar=true);

        //- Construct from IOobject and components.
        //  Boundary is added using addPatches() member function
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& owner,
            labelList&& neighbour,
            const bool syncPar = true
        );

        //- Construct without boundary with cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );

        //- Construct from cell shapes
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const wordList& boundaryPatchTypes,
            const word& defaultBoundaryPatchName,
            const word& defaultBoundaryPatchType,
            const wordList& boundaryPatchPhysicalTypes,
            const bool syncPar = true
        );

        //- Construct from cell shapes, with patch information in dictionary
        //- format.
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const PtrList<dictionary>& boundaryDicts,
            const word& defaultBoundaryPatchName,
            const word& defaultBoundaryPatchType,
            const bool syncPar = true
        );


    //- Destructor
    virtual ~polyMesh();


    // Member Functions

        // Database

            //- Override the objectRegistry dbDir for a single-region case
            virtual const fileName& dbDir() const;

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return the current instance directory for points
            //  Used in the construction of geometric mesh data dependent
            //  on points
            const fileName& pointsInstance() const;

            //- Return the current instance directory for faces
            const fileName& facesInstance() const;

            //- Set the instance for mesh files
            void setInstance
            (
                const fileName& instance,
                const IOobjectOption::writeOption wOpt = IOobject::AUTO_WRITE
            );


        // Regions

            //- Return the local mesh directory name (eg, "polyMesh")
            //- after applying filter for defaultRegion
            static fileName meshDir(const word& region);

            //- The mesh region name or word::null if polyMesh::defaultRegion
            static const word& regionName(const word& region);

            //- The mesh region name or word::null if polyMesh::defaultRegion
            const word& regionName() const;


        // Access

            //- Const reference to the mesh and solver state data
            virtual const meshState& data() const noexcept { return data_; }

            //- Reference to the mesh and solver state data
            virtual meshState& data() noexcept { return data_; }

            //- Return raw points
            virtual const pointField& points() const;

            //- Return true if io is up-to-date with points
            virtual bool upToDatePoints(const regIOobject& io) const;

            //- Set io to be up-to-date with points
            virtual void setUpToDatePoints(regIOobject& io) const;

            //- Return raw faces
            virtual const faceList& faces() const;

            //- Return face owner
            virtual const labelList& faceOwner() const;

            //- Return face neighbour
            virtual const labelList& faceNeighbour() const;

            //- Return old points (mesh motion)
            virtual const pointField& oldPoints() const;

            //- Return old cellCentres (mesh motion)
            virtual const pointField& oldCellCentres() const;

            //- Return boundary mesh
            const polyBoundaryMesh& boundaryMesh() const noexcept
            {
                return boundary_;
            }

            //- Return mesh bounding box
            const boundBox& bounds() const noexcept
            {
                return bounds_;
            }

            //- Return the vector of geometric directions in mesh.
            //  Defined according to the presence of empty and wedge patches.
            //  1 indicates unconstrained direction and -1 a constrained
            //  direction.
            const Vector<label>& geometricD() const;

            //- Return the number of valid geometric dimensions in the mesh
            label nGeometricD() const;

            //- Return the vector of solved-for directions in mesh.
            //  Differs from geometricD in that it includes for wedge cases
            //  the circumferential direction in case of swirl.
            //  1 indicates valid direction and -1 an invalid direction.
            const Vector<label>& solutionD() const;

            //- Return the number of valid solved-for dimensions in the mesh
            label nSolutionD() const;

            //- Return the tetBasePtIs
            const labelIOList& tetBasePtIs() const;

            //- Return the cell search tree
            const indexedOctree<treeDataCell>& cellTree() const;

            //- Return point zone mesh
            const pointZoneMesh& pointZones() const noexcept
            {
                return pointZones_;
            }

            //- Return face zone mesh
            const faceZoneMesh& faceZones() const noexcept
            {
                return faceZones_;
            }

            //- Return cell zone mesh
            const cellZoneMesh& cellZones() const noexcept
            {
                return cellZones_;
            }


        // Database

            //- Return the object registry
            const objectRegistry& thisDb() const noexcept
            {
                return *this;
            }


        // Parallel

            //- The communicator used for parallel communication
            label comm() const noexcept { return comm_; }

            //- The communicator used for parallel communication
            label& comm() noexcept { return comm_; }

            //- Is demand-driven parallel info available?
            bool hasGlobalData() const noexcept;

            //- Return parallel info (demand-driven)
            const globalMeshData& globalData() const;


        // Mesh motion

            //- Is mesh dynamic
            virtual bool dynamic() const
            {
                return false;
            }

            //- Is mesh moving
            bool moving() const noexcept
            {
                return moving_;
            }

            //- Set the mesh to be moving
            bool moving(const bool on) noexcept
            {
                bool old(moving_);
                moving_ = on;
                return old;
            }

            //- Is mesh topology changing
            bool topoChanging() const noexcept
            {
                return topoChanging_;
            }

            //- Set the mesh topology to be changing
            bool topoChanging(const bool on) noexcept
            {
                bool old(topoChanging_);
                topoChanging_ = on;
                return old;
            }

            //- Is mesh changing (topology changing and/or moving)
            bool changing() const noexcept
            {
                return (moving() || topoChanging());
            }

            //- Move points
            virtual void movePoints(const pointField&);

            //- Reset motion
            void resetMotion() const;


        // Topological change

            //- Return non-const access to the pointZones
            pointZoneMesh& pointZones() noexcept
            {
                return pointZones_;
            }

            //- Return non-const access to the faceZones
            faceZoneMesh& faceZones() noexcept
            {
                return faceZones_;
            }

            //- Return non-const access to the cellZones
            cellZoneMesh& cellZones() noexcept
            {
                return cellZones_;
            }

            //- Add boundary patches
            void addPatches
            (
                polyPatchList& plist,
                const bool validBoundary = true
            );

            //- Add boundary patches
            void addPatches
            (
                const List<polyPatch*>& p,
                const bool validBoundary = true
            );

            //- Add mesh zones
            void addZones
            (
                PtrList<pointZone>&& pz,
                PtrList<faceZone>&& fz,
                PtrList<cellZone>&& cz
            );

            //- Add mesh zones
            void addZones
            (
                const List<pointZone*>& pz,
                const List<faceZone*>& fz,
                const List<cellZone*>& cz
            );

            //- Initialise all non-demand-driven data
            virtual bool init(const bool doInit);

            //- Update the mesh based on the mesh files saved in
            //  time directories
            virtual readUpdateState readUpdate();

            //- Update the mesh corresponding to given map
            virtual void updateMesh(const mapPolyMesh& mpm);

            //- Remove boundary patches
            void removeBoundary();

            //- Reset mesh primitive data. Assumes all patch info correct
            //  (so does e.g. parallel communication). If not use
            //  validBoundary=false
            //
            // \note Only autoPtr parameters that test as good() are used
            //    for resetting, otherwise the existing entries are left
            //    untouched.
            void resetPrimitives
            (
                autoPtr<pointField>&& points,
                autoPtr<faceList>&& faces,
                autoPtr<labelList>&& owner,
                autoPtr<labelList>&& neighbour,
                const labelUList& patchSizes,
                const labelUList& patchStarts,
                const bool validBoundary = true
            );


        // Storage management

            //- Clear geometry
            void clearGeom();

            //- Update geometry points; keep topology.
            //- Optionally with new face decomposition
            void updateGeomPoints
            (
                pointIOField&& newPoints,
                autoPtr<labelIOList>& newTetBasePtIsPtr
            );

            //- Clear addressing
            void clearAddressing(const bool isMeshUpdate = false);

            //- Clear all geometry and addressing
            void clearOut(const bool isMeshUpdate = false);

            //- Clear primitive data (points, faces and cells)
            void clearPrimitives();

            //- Clear tet base points
            void clearTetBasePtIs();

            //- Clear cell tree data
            void clearCellTree();

            //- Remove all files from mesh instance
            void removeFiles(const fileName& instanceDir) const;

            //- Remove all files from mesh instance()
            void removeFiles() const;


            bool hasTetBasePtIs() const { return bool(tetBasePtIsPtr_); }


        // Geometric checks. Selectively override primitiveMesh functionality.

            //- Check non-orthogonality
            virtual bool checkFaceOrthogonality
            (
                const bool report = false,
                labelHashSet* setPtr = nullptr
            ) const;

            //- Check face skewness
            virtual bool checkFaceSkewness
            (
                const bool report = false,
                labelHashSet* setPtr = nullptr
            ) const;

            //- Check edge alignment for 1D/2D cases
            virtual bool checkEdgeAlignment
            (
                const bool report,
                const Vector<label>& directions,
                labelHashSet* setPtr
            ) const;

            virtual bool checkCellDeterminant
            (
                const bool report,
                labelHashSet* setPtr
            ) const;

            //- Check mesh motion for correctness given motion points
            virtual bool checkMeshMotion
            (
                const pointField& newPoints,
                const bool report = false,
                const bool detailedReport = false
            ) const;

            //- Check for face weights
            virtual bool checkFaceWeight
            (
                const bool report,
                const scalar minWeight = 0.05,
                labelHashSet* setPtr = nullptr
            ) const;

            //- Check for neighbouring cell volumes
            virtual bool checkVolRatio
            (
                const bool report,
                const scalar minRatio = 0.01,
                labelHashSet* setPtr = nullptr
            ) const;


        // Position search functions

            //- Find the cell, tetFacei and tetPti for point p
            void findCellFacePt
            (
                const point& p,
                label& celli,
                label& tetFacei,
                label& tetPti
            ) const;

            //- Find the tetFacei and tetPti for point p in celli.
            //  tetFacei and tetPtI are set to -1 if not found
            void findTetFacePt
            (
                const label celli,
                const point& p,
                label& tetFacei,
                label& tetPti
            ) const;

            //- Test if point p is in the celli
            bool pointInCell
            (
                const point& p,
                label celli,
                const cellDecomposition = CELL_TETS
            ) const;

            //- Find cell enclosing this location and return index
            //  If not found -1 is returned
            label findCell
            (
                const point& p,
                const cellDecomposition = CELL_TETS
            ) const;


    // Write

        //- Write items held in the objectRegistry. Normally includes mesh
        //- components (points, faces, etc) and any registered fields.
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc = true
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
