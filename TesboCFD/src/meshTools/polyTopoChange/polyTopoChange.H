/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022,2024 OpenCFD Ltd.
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
    Foam::polyTopoChange

Description
    Direct mesh changes based on v1.3 polyTopoChange syntax.

    Instead of recording changes and executing them all in one go (as did
    v1.3 polyTopoChange) this class actually holds the current
    points/faces/cells and does the change immediately.
    It can be asked to compress out all unused points/faces/cells and
    renumber everything to be consistent.

    Note:
    - polyTopoChange can be copied.
    - adding a face using non-existing cells causes all intermediate cells
    to be added. So always first add cells/points and then faces.
    (or set strict checking)
    - strict checking:
        - any added/modified face can only use already existing vertices
        - any added face can only use already existing cells
        - no item can be removed more than once.
    - removed cell: cell set to 0 faces.
    - removed face: face set to 0 vertices.
    - removed point: coordinate set to vector::max (VGREAT,VGREAT,VGREAT).
    Note that this might give problems if this value is used already.
    To see if point is equal to above value we don't use == (which might give
    problems with roundoff error) but instead compare the individual component
    with >.
    - coupled patches: the reorderCoupledFaces routine (borrowed from
    the couplePatches utility) reorders coupled patch faces and
    uses the cyclicPolyPatch,processorPolyPatch functionality.
    - zones are assumed to be non-overlapping by default. If desired to be
    overlapping either set the multiZone to true when calling
    modifyCell|Face|Point or use the variants of addCell|Face|Point
    and modifyCell|Face|Point that take a list of zones.
    - if overlapping zones:
        - 'main' zone is the lowest numbered zone. -1 means no zones.
        - 'additional' zones are stored in incremental ordering (and cannot
          contain -1)

SourceFiles
    polyTopoChange.C
    polyTopoChangeI.H
    polyTopoChangeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef polyTopoChange_H
#define polyTopoChange_H

#include "DynamicList.H"
#include "labelList.H"
#include "pointField.H"
#include "Map.H"
#include "HashSet.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class face;
class primitiveMesh;
class polyMesh;
class Time;
class fileName;
class polyBoundaryMesh;
class polyPatch;
class dictionary;
class topoAction;
class objectMap;
class IOobject;
class mapPolyMesh;
template<class T> class CompactListList;

/*---------------------------------------------------------------------------*\
                       Class polyTopoChange Declaration
\*---------------------------------------------------------------------------*/

class polyTopoChange
{
    // Private data

        //- Whether to allow referencing illegal points/cells/faces
        //  when adding/removing data.
        bool strict_;


        // Patches

            //- Number of patches
            label nPatches_;


        // Points

            //- Current point set
            DynamicList<point> points_;

            //- Original point label (or masterpoint for added points)
            DynamicList<label> pointMap_;

            //- For all original and added points contains new point label.
            //- (used to map return value of addPoint to new mesh point)
            DynamicList<label> reversePointMap_;

            //- Zone of point
            Map<label> pointZone_;

            //- Additional zones of points
            DynamicList<labelList> pointAdditionalZones_;

            //- Retired points
            labelHashSet retiredPoints_;


        // Faces

            //- Current faceList
            DynamicList<face> faces_;

            //- Patch for every external face (-1 for internal faces)
            DynamicList<label> region_;

            //- Owner for all faces
            DynamicList<label> faceOwner_;

            //- Neighbour for internal faces (-1 for external faces)
            DynamicList<label> faceNeighbour_;

            //- Original face label. Or master face for added-from-faces;
            //  -1 for faces added-from-edge or added-from-point)
            DynamicList<label> faceMap_;

            //- For all original and added faces contains new face label
            //  (used to map return value of addFace to new mesh face)
            DynamicList<label> reverseFaceMap_;

            //- Faces added from point (corresponding faceMap_ will
            //  be -1)
            Map<label> faceFromPoint_;

            //- Faces added from edge (corresponding faceMap_ will
            //  be -1)
            Map<label> faceFromEdge_;

            //- In mapping whether to reverse the flux.
            bitSet flipFaceFlux_;

            //- Zone of face
            Map<label> faceZone_;

            //- Orientation of face in zone
            bitSet faceZoneFlip_;

            //- Additional zones of face (stored as signed face)
            DynamicList<labelList> faceAdditionalZones_;

            //- Active faces
            label nActiveFaces_;


        // Cells

            //- Original cell label or master cell for added-from-cell;
            //  -1 for cells added from face or edge.
            DynamicList<label> cellMap_;

            //- For all original and added cells contains new cell label
            //  (used to map return value of addCell to new mesh cell)
            DynamicList<label> reverseCellMap_;

            //- Cells added from point
            Map<label> cellFromPoint_;

            //- Cells added from edge
            Map<label> cellFromEdge_;

            //- Cells added from face
            Map<label> cellFromFace_;

            //- First zone of cell (or -1)
            DynamicList<label> cellZone_;

            //- Additional zones of cell
            DynamicList<labelList> cellAdditionalZones_;


    // Private Member Functions

        //- Reorder contents of container according to oldToNew map
        template<class Type>
        static void reorder
        (
            const labelUList& oldToNew,
            DynamicList<Type>& lst
        );

        template<class Type>
        static void renumberKey
        (
            const labelUList& oldToNew,
            Map<Type>& map
        );

        //- Renumber elements of container according to oldToNew map
        static void renumber
        (
            const labelUList& oldToNew,
            labelHashSet& labels
        );

        //- Special handling of reverse maps which have <-1 in them
        static void renumberReverseMap
        (
            const labelUList& oldToNew,
            DynamicList<label>& elems
        );

        //- Renumber & compact elements of list according to map
        static void renumberCompact
        (
            const labelUList& oldToNew,
            labelList& elems
        );

        //- Count number of added and removed quantities from maps.
        static void countMap
        (
            const labelUList& map,
            const labelUList& reverseMap,
            label& nAdd,
            label& nInflate,
            label& nMerge,
            label& nRemove
        );

        //- Print some stats about mesh
        static void writeMeshStats(const polyMesh& mesh, Ostream& os);

        //- Calculate object maps. Requires reverseMap to have destination
        //  to be marked with <-1.
        static void getMergeSets
        (
            const labelUList& reverseCellMap,
            const labelUList& cellMap,
            List<objectMap>& cellsFromCells
        );

        //- Are all face vertices valid
        bool hasValidPoints(const face& f) const;

        //- Return face points
        pointField facePoints(const face& f) const;

        //- Check inputs to modFace or addFace
        void checkFace
        (
            const face& f,
            const label facei,
            const label own,
            const label nei,
            const label patchi,
            const label zoneI
        ) const;

        //- Construct cells (in packed storage)
        void makeCells
        (
            const label nActiveFaces,
            labelList& cellFaces,
            labelList& cellFaceOffsets
        ) const;

        //- Construct cellCells (in packed storage)
        void makeCellCells
        (
            const label nActiveFaces,
            CompactListList<label>& cellCells
        ) const;

        //- Cell ordering (bandCompression). Returns number of remaining cells.
        label getCellOrder
        (
            const CompactListList<label>& cellCellAddressing,
            labelList& oldToNew
        ) const;

        //- Do upper-triangular ordering and patch ordering.
        void getFaceOrder
        (
            const label nActiveFaces,
            const labelUList& cellFaces,
            const labelUList& cellFaceOffsets,

            labelList& oldToNew,
            labelList& patchSizes,
            labelList& patchStarts
        ) const;

        //- Compact and reorder faces according to map
        void reorderCompactFaces
        (
            const label newSize,
            const labelUList& oldToNew
        );

        //- Remove all unused/removed points/faces/cells and update
        //- face ordering (always), cell ordering (bandcompression,
        //- orderCells=true),
        //- point ordering (sorted into internal and boundary points,
        //- orderPoints=true)
        void compact
        (
            const bool orderCells,
            const bool orderPoints,
            label& nInternalPoints,
            labelList& patchSizes,
            labelList& patchStarts
        );

        //- Select either internal or external faces out of faceLabels
        static labelList selectFaces
        (
            const primitiveMesh& mesh,
            const labelUList& faceLabels,
            const bool internalFacesOnly
        );

        //- Calculate mapping for patchpoints only
        void calcPatchPointMap
        (
            const UList<Map<label>>& oldPatchMeshPointMaps,
            const labelUList& patchMap,
            const polyBoundaryMesh& boundary,
            labelListList& patchPointMap
        ) const;

        void calcFaceInflationMaps
        (
            const polyMesh& mesh,
            List<objectMap>& facesFromPoints,
            List<objectMap>& facesFromEdges,
            List<objectMap>& facesFromFaces
        ) const;

        void calcCellInflationMaps
        (
            const polyMesh& mesh,
            List<objectMap>& cellsFromPoints,
            List<objectMap>& cellsFromEdges,
            List<objectMap>& cellsFromFaces,
            List<objectMap>& cellsFromCells
        ) const;

        void resetZones
        (
            const polyMesh& mesh,   // mesh to get existing info from
            polyMesh& newMesh,      // mesh to change zones on
            labelListList& pointZoneMap,
            labelListList& faceZoneFaceMap,
            labelListList& cellZoneMap
        ) const;

        void calcFaceZonePointMap
        (
            const polyMesh& mesh,
            const UList<Map<label>>& oldFaceZoneMeshPointMaps,
            labelListList& faceZonePointMap
        ) const;


        // Coupling

            //- Do all coupled patch face reordering
            void reorderCoupledFaces
            (
                const bool syncParallel,
                const polyBoundaryMesh& boundary,
                const labelUList& patchMap,     // new back to old patch labels
                const labelUList& patchStarts,
                const labelUList& patchSizes,
                const pointField& points
            );

        void compactAndReorder
        (
            const polyMesh& mesh,
            const labelUList& patchMap,      // from new to old patch
            const bool syncParallel,
            const bool orderCells,
            const bool orderPoints,
            label& nInternalPoints,
            pointField& newPoints,
            labelList& patchSizes,
            labelList& patchStarts,
            List<objectMap>& pointsFromPoints,
            List<objectMap>& facesFromPoints,
            List<objectMap>& facesFromEdges,
            List<objectMap>& facesFromFaces,
            List<objectMap>& cellsFromPoints,
            List<objectMap>& cellsFromEdges,
            List<objectMap>& cellsFromFaces,
            List<objectMap>& cellsFromCells,
            List<Map<label>>& oldPatchMeshPointMaps,
            labelList& oldPatchNMeshPoints,
            labelList& oldPatchStarts,
            List<Map<label>>& oldFaceZoneMeshPointMaps
        );

public:

    //- Runtime type information
    ClassName("polyTopoChange");



    // Constructors

        //- Construct without mesh. Either specify nPatches or use
        //- setNumPatches before trying to make a mesh (makeMesh, changeMesh)
        polyTopoChange(const label nPatches, const bool strict = true);

        //- Construct from mesh. Adds all points/face/cells from mesh.
        polyTopoChange(const polyMesh& mesh, const bool strict = true);


    // Member Functions

        // Access

            //- Points. Shrunk after constructing mesh (or calling of compact())
            const DynamicList<point>& points() const
            {
                return points_;
            }

            const DynamicList<face>& faces() const
            {
                return faces_;
            }

            const DynamicList<label>& region() const
            {
                return region_;
            }

            const DynamicList<label>& faceOwner() const
            {
                return faceOwner_;
            }

            const DynamicList<label>& faceNeighbour() const
            {
                return faceNeighbour_;
            }

            //- Is point removed?
            //- Considered removed if point is GREAT.
            inline bool pointRemoved(const label pointi) const;

            //- Is face removed?
            //- Considered removed if face is empty
            inline bool faceRemoved(const label facei) const;

            //- Is cell removed?
            //- Considered removed if the cellMap is -2
            inline bool cellRemoved(const label celli) const;


        // Edit

            //- Clear all storage
            void clear();

            //- Add all points/faces/cells of mesh. Additional offset for patch
            //- or zone ids.
            void addMesh
            (
                const polyMesh& mesh,
                const labelUList& patchMap,
                const labelUList& pointZoneMap,
                const labelUList& faceZoneMap,
                const labelUList& cellZoneMap
            );

            //- Explicitly pre-size the dynamic storage for expected mesh
            //- size for if construct-without-mesh
            void setCapacity
            (
                const label nPoints,
                const label nFaces,
                const label nCells
            );

            //- Shrink storage (does not remove any elements; just compacts
            //- dynamic lists
            void shrink();

            //- Move all points. Incompatible with other topology changes.
            void movePoints(const pointField& newPoints);

            //- For compatibility with polyTopoChange: set topological action.
            label setAction(const topoAction& action);

            //- Add point. Return new point label.
            //  Notes:
            //  - masterPointID can be < 0 (appended points)
            //  - inCell = false: add retired point (to end of point list)
            label addPoint
            (
                const point& pt,
                const label masterPointID,
                const label zoneID,
                const bool inCell
            );

            //- Add point. Return new point label.
            //  Notes:
            //  - masterPointID can be < 0 (appended points)
            //  - inCell = false: add retired point (to end of point list)
            label addPoint
            (
                const point& pt,
                const label masterPointID,
                const labelUList& zoneIDs,
                const bool inCell
            );

            //- Modify coordinate.
            //  Notes:
            //  - zoneID = +ve (add to zoneID), -ve (remove from zones)
            //  - inCell = false: add retired point (to end of point list)
            void modifyPoint
            (
                const label pointi,
                const point& pt,
                const label zoneID,
                const bool inCell,
                const bool multiZone = false
            );

            //- Modify coordinate.
            //  Notes:
            //  - zoneIDs = set pointZones
            //  - inCell = false: add retired point (to end of point list)
            void modifyPoint
            (
                const label pointi,
                const point& pt,
                const labelUList& zoneIDs,
                const bool inCell
            );

            //- Remove/merge point.
            void removePoint(const label pointi, const label mergePointi);

            //- Get current cellZone(s). Return number of zones.
            label pointZones
            (
                const label pointi,
                DynamicList<label>& zones
            ) const;

            //- Add face to cells. Return new face label.
            //  own,nei<0, zoneID>=0 : add inactive face (to end of face list)
            label addFace
            (
                const face& f,
                const label own,
                const label nei,
                const label masterPointID,
                const label masterEdgeID,
                const label masterFaceID,
                const bool flipFaceFlux,
                const label patchID,
                const label zoneID,
                const bool zoneFlip
            );

            //- Add face to cells. Return new face label.
            //  own,nei<0, zoneID>=0 : add inactive face (to end of face list)
            label addFace
            (
                const face& f,
                const label own,
                const label nei,
                const label masterPointID,
                const label masterEdgeID,
                const label masterFaceID,
                const bool flipFaceFlux,
                const label patchID,
                const labelUList& zoneIDs,
                const UList<bool>& zoneFlips
            );

            //- Modify vertices or cell of face.
            void modifyFace
            (
                const face& f,
                const label facei,
                const label own,
                const label nei,
                const bool flipFaceFlux,
                const label patchID,
                const label zoneID,
                const bool zoneFlip,
                const bool multiZone = false
            );

            //- Modify vertices or cell of face.
            void modifyFace
            (
                const face& f,
                const label facei,
                const label own,
                const label nei,
                const bool flipFaceFlux,
                const label patchID,
                const labelUList& zoneIDs,
                const UList<bool>& zoneFlips
            );

            //- Remove/merge face.
            void removeFace(const label facei, const label mergeFacei);

            //- Get current faceZone(s). Return number of zones.
            label faceZones
            (
                const label facei,
                DynamicList<label>& zones,
                DynamicList<bool>& flips
            ) const;

            //- Add cell. Return new cell label.
            label addCell
            (
                const label masterPointID,
                const label masterEdgeID,
                const label masterFaceID,
                const label masterCellID,
                const label zoneID
            );

            //- Add zoned cell (zones cannot be -1). Return new cell label.
            label addCell
            (
                const label masterPointID,
                const label masterEdgeID,
                const label masterFaceID,
                const label masterCellID,
                const labelUList& zoneIDs
            );

            //- Modify zone of cell. Optionally add to zone.
            void modifyCell
            (
                const label celli,
                const label zoneID,
                const bool multiZone = false
            );

            //- Set zones of cell
            void modifyCell(const label celli, const labelUList& zoneIDs);

            //- Remove/merge cell.
            void removeCell(const label celli, const label mergeCelli);

            //- Get current cellZone(s). Return number of zones.
            label cellZones(const label celli, DynamicList<label>& zones) const;


            //- Explicitly set the number of patches if construct-without-mesh
            //- used.
            inline void setNumPatches(const label nPatches);

        // Other

            //- Inplace changes mesh without change of patches.
            //  Adapts patch start/end and by default does parallel matching.
            //  Clears all data. Returns map.
            //  inflate = true : keep old mesh points. Put new points into the
            //  returned map (preMotionPoints) so we can use inflation. Any
            //  points out of nothing (appended points) are vector::zero.
            //  inflate = false: set mesh points directly. Empty preMotionPoints
            //  in the map.
            //  orderCells :  whether to order the cells (see bandCompression.H)
            //  orderPoints : whether to order the points into internal first
            //  followed by boundary points. This is not fully consistent
            //  with upper-triangular ordering of points and edges so
            //  is only done when explicitly asked for.
            autoPtr<mapPolyMesh> changeMesh
            (
                polyMesh& mesh,
                const labelUList& patchMap,
                const bool inflate,
                const bool syncParallel = true,
                const bool orderCells = false,
                const bool orderPoints = false
            );

            //- Without patch mapping
            autoPtr<mapPolyMesh> changeMesh
            (
                polyMesh& mesh,
                const bool inflate,
                const bool syncParallel = true,
                const bool orderCells = false,
                const bool orderPoints = false
            );

            //- Create new mesh with old mesh patches. Additional dictionaries
            //  (fv* etc) read according to IO flags
            template<class Type>
            autoPtr<mapPolyMesh> makeMesh
            (
                autoPtr<Type>& newMesh,
                const IOobject& io,
                const polyMesh& mesh,
                const labelUList& patchMap,
                const bool syncParallel = true,
                const bool orderCells = false,
                const bool orderPoints = false
            );

            //- Create new mesh with old mesh patches. Additional dictionaries
            //  (fv* etc) read according to IO flags
            template<class Type>
            autoPtr<mapPolyMesh> makeMesh
            (
                autoPtr<Type>& newMesh,
                const IOobject& io,
                const polyMesh& mesh,
                const bool syncParallel = true,
                const bool orderCells = false,
                const bool orderPoints = false
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyTopoChangeI.H"

#ifdef NoRepository
    #include "polyTopoChangeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
