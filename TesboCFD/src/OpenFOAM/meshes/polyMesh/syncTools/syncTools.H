/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::syncTools

Description
    Various tools to aid synchronizing lists across coupled patches. WIP.

    Require
    - combineOperator (e.g. sumEqOp - not sumOp!) that is defined for the
      type be defined.
    - null value which gets overridden by any valid value.
    - transform function

SourceFiles
    syncTools.C
    syncToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_syncTools_H
#define Foam_syncTools_H

#include "Pstream.H"
#include "edgeHashes.H"
#include "bitSet.H"
#include "polyMesh.H"
#include "coupledPolyPatch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class syncTools Declaration
\*---------------------------------------------------------------------------*/

class syncTools
{
    // Private Member Functions

        //- Combine val with existing value in pointValues map at given index
        //  No communication
        template<class T, class CombineOp>
        static void combine
        (
            Map<T>& pointValues,
            const CombineOp& cop,
            const label index,
            const T& val
        );

        //- Combine val with existing value in edgeValues at edge index
        //  No communication
        template<class T, class CombineOp>
        static void combine
        (
            EdgeMap<T>& edgeValues,
            const CombineOp& cop,
            const edge& index,
            const T& val
        );


public:

        // Basic routines with user-supplied transformation.
        // Preferably use specialisations below.

            //- Synchronize values on selected points.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointMap
            (
                const polyMesh& mesh,
                Map<T>& pointValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on selected edges.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp, class TransformOp>
            static void syncEdgeMap
            (
                const polyMesh& mesh,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on all mesh points.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on selected mesh points.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                const labelUList& meshPoints,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on all mesh edges.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp, class TransformOp, class FlipOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top,
                const FlipOp& fop
            );

            //- Synchronize values on selected mesh edges.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp, class TransformOp, class FlipOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                const labelUList& meshEdges,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top,
                const FlipOp& fop
            );

            //- Synchronize values on boundary faces only.
            template<class T, class CombineOp, class TransformOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                const CombineOp& cop,
                const TransformOp& top,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            );


        // Synchronise point-wise data

            //- Synchronize values on all mesh points.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    pointValues,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on all mesh points.
            //  Communication if UPstream::parRun() == true.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on selected mesh points.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                const labelUList& meshPoints,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    pointValues,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected mesh points.
            //  Communication if UPstream::parRun() == true.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                const labelUList& meshPoints,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }


        // Synchronise edge-wise data

            //- Synchronize values on all mesh edges.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    edgeValues,
                    cop,
                    nullValue,
                    mapDistribute::transform(),
                    identityOp()  // No flipping
                );
            }

            //- Synchronize locations on all mesh edges.
            //  Communication if UPstream::parRun() == true.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition(),
                    identityOp()  // No flipping
                );
            }

            //- Synchronize values on selected mesh edges.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                const labelUList& meshEdges,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    meshEdges,
                    edgeValues,
                    cop,
                    nullValue,
                    mapDistribute::transform(),
                    identityOp()  // No flipping
                );
            }

            //- Synchronize locations on selected mesh edges.
            //  Communication if UPstream::parRun() == true.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                const labelUList& meshEdges,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    meshEdges,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition(),
                    identityOp()  // No flipping
                );
            }



        // Synchronise face-wise data

            //- Synchronize values on boundary faces only.
            template<class T, class CombineOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    faceValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on boundary faces only.
            template<class CombineOp>
            static void syncBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<point>& positions,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    positions,
                    cop,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on all mesh faces.
            template<class T, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                const CombineOp& cop,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            )
            {
                SubList<T> bndValues
                (
                    faceValues,
                    mesh.nBoundaryFaces(),
                    mesh.nInternalFaces()
                );

                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    cop,
                    mapDistribute::transform(),
                    parRun
                );
            }

            //- Synchronize locations on all mesh faces.
            template<class CombineOp>
            static void syncFacePositions
            (
                const polyMesh& mesh,
                UList<point>& positions,
                const CombineOp& cop,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            )
            {
                SubList<point> bndValues
                (
                    positions,
                    mesh.nBoundaryFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    cop,
                    mapDistribute::transformPosition(),
                    parRun
                );
            }

            //- Swap coupled boundary face values. Uses eqOp
            template<class T>
            static void swapBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    faceValues,
                    eqOp<T>(),
                    mapDistribute::transform(),
                    parRun
                );
            }

            //- Swap coupled positions. Uses eqOp
            static void swapBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<point>& positions,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    positions,
                    eqOp<point>(),
                    mapDistribute::transformPosition(),
                    parRun
                );
            }

            //- Swap coupled face values. Uses eqOp
            template<class T>
            static void swapFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            )
            {
                SubList<T> bndValues
                (
                    faceValues,
                    mesh.nBoundaryFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    eqOp<T>(),
                    mapDistribute::transform(),
                    parRun
                );
            }

            //- Extract and swap to obtain neighbour cell values
            //- for all boundary faces
            template<class T>
            static void swapBoundaryCellList
            (
                const polyMesh& mesh,
                const UList<T>& cellData,
                List<T>& neighbourCellData,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            );

            //- Extract and swap to obtain neighbour cell positions
            //- for all boundary faces
            static void swapBoundaryCellPositions
            (
                const polyMesh& mesh,
                const UList<point>& cellData,
                List<point>& neighbourCellData,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            );


            //- Return neighbour cell values for all boundary faces
            //- by swapping via boundary faces
            template<class T>
            FOAM_NODISCARD static List<T> swapBoundaryCellList
            (
                const polyMesh& mesh,
                const UList<T>& cellData,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            )
            {
                List<T> nbrCellData;
                swapBoundaryCellList(mesh, cellData, nbrCellData, parRun);
                return nbrCellData;
            }

            //- Return neighbour cell positions for all boundary faces
            //- by swapping via boundary faces
            FOAM_NODISCARD static List<point> swapBoundaryCellPositions
            (
                const polyMesh& mesh,
                const UList<point>& cellData,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            )
            {
                List<point> nbrCellData;
                swapBoundaryCellPositions(mesh, cellData, nbrCellData, parRun);
                return nbrCellData;
            }


        // Sparse versions

            //- Synchronize values on selected points.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp>
            static void syncPointMap
            (
                const polyMesh& mesh,
                Map<T>& pointValues,
                const CombineOp& cop
            )
            {
                syncPointMap
                (
                    mesh,
                    pointValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected points.
            //  Communication if UPstream::parRun() == true.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                Map<point>& positions,
                const CombineOp& cop
            )
            {
                syncPointMap
                (
                    mesh,
                    positions,
                    cop,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on selected edges.
            //- Edges are represented by the two vertices that make it up
            //- so global edges never get constructed.
            //  Communication if UPstream::parRun() == true.
            template<class T, class CombineOp>
            static void syncEdgeMap
            (
                const polyMesh& mesh,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop
            )
            {
                syncEdgeMap
                (
                    mesh,
                    edgeValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected edges.
            //  Communication if UPstream::parRun() == true.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                EdgeMap<point>& edgePositions,
                const CombineOp& cop
            )
            {
                syncEdgeMap
                (
                    mesh,
                    edgePositions,
                    cop,
                    mapDistribute::transformPosition()
                );
            }


        // PackedList versions

            //- Synchronize face values from PackedList/bitSet
            //
            //  \param mesh The mesh
            //  \param isBoundaryOnly True if faceValues are for the boundary
            //      only and not the entire mesh. This determines the face
            //      offset when accessing values.
            //  \param faceValues The face values to synchronize
            //  \param cop The combine operation
            //  \param parRun Allow parallel communication
            template<unsigned Width, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                const bool isBoundaryOnly,
                PackedList<Width>& faceValues,
                const CombineOp& cop,
                const bool parRun = UPstream::parRun()
            );

            //- Synchronize mesh face values from PackedList/bitSet
            template<unsigned Width, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues,
                const CombineOp& cop,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            );

            //- Synchronize boundary face values from PackedList/bitSet
            template<unsigned Width, class CombineOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues,
                const CombineOp& cop,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            );

            //- Swap coupled face values. Uses eqOp
            template<unsigned Width>
            static void swapFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            );

            //- Swap coupled boundary face values. Uses eqOp
            template<unsigned Width>
            static void swapBoundaryFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues,
                //! Allow parallel communication
                const bool parRun = UPstream::parRun()
            );

            template<unsigned Width, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                PackedList<Width>& pointValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

            template<unsigned Width, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                PackedList<Width>& edgeValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );


        // Other

            //- Get per point whether it is uncoupled or a master of a
            //- coupled set of points
            static bitSet getMasterPoints(const polyMesh& mesh);

            //- Get per edge whether it is uncoupled or a master of a
            //- coupled set of edges
            static bitSet getMasterEdges(const polyMesh& mesh);

            //- Get per face whether it is uncoupled or a master of a
            //- coupled set of faces
            static bitSet getMasterFaces(const polyMesh& mesh);

            //- Get per face whether it is internal or a master of a
            //- coupled set of faces
            static bitSet getInternalOrMasterFaces(const polyMesh& mesh);

            //- Get per face whether it is internal or coupled
            static bitSet getInternalOrCoupledFaces(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "syncToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
