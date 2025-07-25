/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020,2024 OpenCFD Ltd.
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
    Foam::addPatchCellLayer

Description
    Adds layers of cells to outside (or inside) of polyMesh. Can optionally
    create stand-alone extruded mesh (addToMesh=false).

    Call setRefinement with offset vector for every patch point and number
    of layers per patch face and number of layers per patch point.
    - offset vector should be zero for any non-manifold point and synchronised
      on coupled points before calling this.
    - offset vector of zero will not add any points.
    - gets supplied the number of extruded layers both per face and per
      point. Usually the point nlayers is the max of surrounding face nlayers.

      point nlayers:
       -  0 : no extrusion. Any surrounding face being extruded becomes 'prism'
       - >0 : should be max of surrounding face nlayers.

    - differing face nlayers: 'termination' : (e.g. from 2 to 4 layers) match
      at original patch face side.

        E.g. 2 boundary faces on patches a,b. 2 layers for a, 3 for b.

    \verbatim
        Was:

               a      b         <- patch of boundary face
            +------+------+
            |      |      |     <- original cells
            +------+------+

        Extrusion:


                            added   added
                            face    cell
                            ----    ----
               a      b                     <- patch of boundary face
            +------+------+  3
            |      |      |         2
            +      +------+  2
            |      |      |         1
            +------+------+  1
            |      |      |         0
            +------+------+  0              <- original boundary faces
            |      |      |                 <- original cells
            +------+------+



        Intrusion:

                            face    cell
                            ----    ----
               a      b                     <- patch of boundary face
            +------+------+  0
            |      |      |         0
            +------+------+  1
            |      |      |         1
            +      +------+  2
            |      |      |         2
            +------+------+  3
            |      |      |                 <- original cells
            +------+------+
    \endverbatim


    - 'side' faces (i.e. on the edge of pp) get the patchID/zoneID of the
    other patch/zone they are connected to (hopefully only 1)


    E.g. 3 boundary faces on patches a,b. b gets extruded, a doesn't.

    \verbatim
           a      b      b        <- patch of boundary face
        +------+------+------+
        |      |      |      |
        |      |      |      |    <- cells
        |      |      |      |
        +------+------+------+


        (shown for extrusion mode only):


               ^      ^           <- wanted extrusion vector (none at far right)
           a   |  b   |  b        <- patch of boundary face
        +------+------+------+
        |      |      |      |    <- cells
        +------+------+------+

                  b
               +------+\ b        1. prism cell added onto second b face since
           a  a|      | ----\          only one side gets extruded.
        +------+------+------+    2. side-face gets patch a, not b.
        |      |      |      |
        +------+------+------+
    \endverbatim


SourceFiles
    addPatchCellLayer.C

\*---------------------------------------------------------------------------*/

#ifndef addPatchCellLayer_H
#define addPatchCellLayer_H

#include "labelList.H"
#include "typeInfo.H"
#include "labelPair.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class polyTopoChange;
class mapPolyMesh;
class primitiveMesh;
class globalIndex;

/*---------------------------------------------------------------------------*\
                      Class addPatchCellLayer Declaration
\*---------------------------------------------------------------------------*/

class addPatchCellLayer
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Add layers to existing mesh or create new mesh
        const bool addToMesh_;

        //- Add layers to outside of mesh or to inside
        const bool extrude_;


        //- For all patchpoints: list of added points (size 0 or nLayers)
        //  First point in list is one nearest to original point in patch,
        //  last one is
        //  - extrude : the new point on the surface
        //  - intrude : the point connecting to the original cell.
        labelListList addedPoints_;

        //- For all patchfaces: list of layer faces.
        //  - empty if no face extruded
        //  - first element is original patch face
        //  - last element is
        //      - extrude : the new boundary face
        //      - intrude : the new internal face to the original cell.
        labelListList layerFaces_;


    // Private Member Functions

        //- Get the face on the other side of the edge.
        static label nbrFace
        (
            const labelListList& edgeFaces,
            const label edgeI,
            const label faceI
        );

        //- Add vertex to face if unique.
        static void addVertex(const label, face&, label& fp);

        bool sameEdgeNeighbour
        (
            const indirectPrimitivePatch& pp,
            const labelListList& globalEdgeFaces,
            const boolList& doneEdge,
            const label thisGlobalFaceI,
            const label nbrGlobalFaceI,
            const label edgeI
        ) const;

        labelPair getEdgeString
        (
            const indirectPrimitivePatch& pp,
            const labelListList& globalEdgeFaces,
            const boolList& doneEdge,
            const label patchFaceI,
            const label globalFaceI
        ) const;


        //- Add face between layer-1 and layer.
        label addSideFace
        (
            const indirectPrimitivePatch&,
            const labelListList& addedCells,

            const face& newFace,
            const label newPatchID,
            const label newZoneI,
            const bool newFlip,
            const label inflateFaceI,

            const label ownFaceI,
            const label nbrFaceI,
            const label meshEdgeI,
            const label layerI,
            const label numEdgeFaces,
            const labelList& meshFaces,
            polyTopoChange&
        ) const;

        //- Extract properties from mesh face
        static void setFaceProps
        (
            const polyMesh&,
            const label,
            label&,
            label&,
            bool&
        );

        //- Extract properties from mesh face in pp edge ordering
        static void setFaceProps
        (
            const polyMesh& mesh,
            const indirectPrimitivePatch& pp,
            const label ppEdgeI,
            const label faceI,

            label& patchI,
            label& zoneI,
            bool& zoneFlip,
            label& inflateFaceI
        );

        //- Mark internal and boundary edges of patch. In mesh edges
        //- since processor might not have pp but does have edge.
        static void markPatchEdges
        (
            const polyMesh& mesh,
            const indirectPrimitivePatch& pp,
            const labelListList& edgeGlobalFaces,
            const labelList& meshEdges,

            bitSet& isPatchEdge,
            bitSet& isPatchBoundaryEdge
        );

        //- For every edge on pp return
        //  - patchEdgeToFace : face (in global indexing) to inflate from
        //  - patchEdgeToPatch : patch (only for boundary edges of pp)
        //  - patchEdgeToZone,flip : zone info
        static void globalEdgeInfo
        (
            const bool zoneFromAnyFace,

            const polyMesh& mesh,
            const globalIndex& globalFaces,
            const labelListList& edgeGlobalFaces,
            const indirectPrimitivePatch& pp,
            const labelList& meshEdges,

            labelList& patchEdgeToFace,  // face (in globalFaces index)
            labelList& patchEdgeToPatch, // patch on face (or -1 for int faces)
            labelList& patchEdgeToZone,  // zone on face
            bitSet& patchEdgeToFlip      // flip orientation on face
        );

        //- No copy construct
        addPatchCellLayer(const addPatchCellLayer&) = delete;

        //- No copy assignment
        void operator=(const addPatchCellLayer&) = delete;


public:

    //- Runtime type information
    ClassName("addPatchCellLayer");


    // Constructors

        //- Construct from mesh.
        explicit addPatchCellLayer
        (
            const polyMesh&,
            const bool addToMesh=true,
            const bool extrude=true
        );


    // Member Functions

        // Access

            //- Added points per patch point.
            const labelListList& addedPoints() const
            {
                return addedPoints_;
            }

            //- Layer faces per patch face. See above.
            const labelListList& layerFaces() const
            {
                return layerFaces_;
            }

            //- Helper: get added cells per patch face.
            //  addedCells[patchFace] is list of cells added.
            //  extrude :
            //      first element : next to original cell
            //      last element  : is the top cell (i.e. the boundary cell)
            //  intrude :
            //      first element : top cell
            //      last element  : next to original cell
            static labelListList addedCells
            (
                const polyMesh&,
                const labelListList& layerFaces
            );

            //- Added cells given current mesh & layerfaces.
            labelListList addedCells() const;


        // Edit

            //- Per patch edge the pp faces (in global indices) using it.
            //  Uses ListOps::uniqueEqOp to remove duplicates. On coupled
            //  faces only selects the one with the correct orientation/flip
            //  (assumes the orientation is opposite on a coupled face pair)
            static labelListList globalEdgeFaces
            (
                const polyMesh&,
                const globalIndex& globalFaces,
                const indirectPrimitivePatch& pp,
                const bitSet& orientation
            );

            //- Per patch edge the pp faces (in global indices) using it.
            //  Uses ListOps::uniqueEqOp to remove duplicates.
            static labelListList globalEdgeFaces
            (
                const polyMesh&,
                const globalIndex& globalFaces,
                const indirectPrimitivePatch& pp
            );

            //- Determine extrude information per patch edge:
            //  - zoneID, zoneFlip :
            //      picks one of the faces that connects to
            //      the edge. For boundary edge only looks
            //      at boundary faces. For internal edge it looks at internal
            //      faces only (zoneFromAnyFace = false) or at any face
            //      (zoneFromAnyFace = true). zoneFlip is consistent with
            //      ordering of pp edge.
            //      Face selected gets stored in inflateFaceID
            //  - patchID :
            //      get patch from any boundary face connected to the
            //      edge. The patch might be a to-be-created processor patch
            //      (patchI >= mesh.boundaryMesh().size()) in which case the
            //      nbrProcToPatch, patchToNbrProc give the correspondence.
            //      nPatches is the new number of patches.
            static void calcExtrudeInfo
            (
                const bool zoneFromAnyFace,

                const polyMesh&,
                const globalIndex& globalFaces,
                const labelListList& globalEdgeFaces,
                const indirectPrimitivePatch& pp,

                labelList& edgePatchID,     // if extruding a patch edge
                label& nPatches,
                Map<label>& nbrProcToPatch,
                Map<label>& patchToNbrProc,
                labelList& edgeZoneID,
                boolList& edgeFlip,
                labelList& inflateFaceID
            );

            //- Play commands into polyTopoChange to create layers on top
            //  of indirectPrimitivePatch (have to be outside faces).
            //  Gets displacement per patch point.
            //  - exposedPatchID : only used if creating a new mesh
            //    (addToMesh=false) gives per pp face the patch the
            //    exposed face should get.
            //  - nPointLayers : number of layers per (patch)point.
            //  - nFaceLayers : number of layers per (patch) face.
            //  - firstDisplacement : displacement per point for first
            //    layer of points (i.e. nearest to original mesh). If zero
            //    do not add point.
            //  Layer thicknesses are calculated to constant geometric
            //  expansion. Use expansionRatio 1 for constant size.
            //  Sets addedPoints_ which is per pp point a list of points
            //  added.
            //  Note: firstDisplacement has to be parallel synchronised before
            //        calling this routine. Only if all procs sharing a point
            //        get a cell should firstDisplacement be <> 0
            //  Note: cells get added from owner cells of patch faces
            //        (instead of e.g. from patch faces)
            void setRefinement
            (
                const globalIndex& globalFaces,
                const labelListList& globalEdgeFaces,
                const scalarField& expansionRatio,
                const indirectPrimitivePatch& pp,
                const bitSet& flip,

                const labelList& sidePatchID,
                const labelList& sideZoneID,
                const boolList& sideFlip,
                const labelList& inflateFaceID,

                const labelList& exposedPatchID,
                const labelList& nFaceLayers,
                const labelList& nPointLayers,
                const vectorField& firstLayerDisp,
                polyTopoChange& meshMod
            );


            //- Add with constant expansion ratio and same nLayers everywhere
            void setRefinement
            (
                const globalIndex& globalFaces,
                const labelListList& globalEdgeFaces,
                const label nLayers,
                const indirectPrimitivePatch& pp,
                const labelList& sidePatchID,
                const vectorField& overallDisplacement,
                polyTopoChange& meshMod
            )
            {
                setRefinement
                (
                    globalFaces,
                    globalEdgeFaces,
                    scalarField(pp.nPoints(), 1.0),     // expansion ration
                    pp,
                    bitSet(pp.size()),                  // flip

                    sidePatchID,
                    labelList(pp.nEdges(), -1),         // zoneID
                    boolList(pp.nEdges(), false),       // zoneFlip
                    labelList(pp.nEdges(), -1),         // inflateFaceID

                    labelList(0),
                    labelList(pp.size(), nLayers),      // nFaceLayers
                    labelList(pp.nPoints(), nLayers),   // nPointLayers
                    overallDisplacement / nLayers,      // firstLayerDisp
                    meshMod
                );
            }


            //- Update any locally stored mesh information. Gets additional
            //  map from new to old patch (since patch needs to be
            //  recreated since has to be on outside).
            void updateMesh
            (
                const mapPolyMesh&,
                const labelList& faceMap,   // new to old patch faces
                const labelList& pointMap   // new to old patch points
            );



            //- Helper: given patch and points on patch that are extruded
            //  (to slave side or master side) find the affected
            //  points. Calculates by walking across faces which vertices on
            //  which face are affected. isDupMeshPoint:
            //      -1  : unaffected
            //      >=0 : should use local duplicate of point
            //            (though it does not tell us whether it should use
            //             slave side or master side)
            static void findDuplicatedPoints
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& pp,
                const bitSet& ppFlip,       // optional orientation on top of pp
                const bitSet& isBlockedFace,// any mesh faces not to be
                                            // traversed.Usually pp.addressing()
                const bitSet& isDupPatchPoint,
                const bool extrude,         // which side to extrude
                faceList& isDupMeshPoint
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
