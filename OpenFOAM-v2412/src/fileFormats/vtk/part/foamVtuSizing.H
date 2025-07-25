/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::vtk::vtuSizing

Description
    Sizing descriptions and routines for transcribing an OpenFOAM volume mesh
    into a VTK unstructured grid, with possible decomposition of polyhedral
    cells into primitive cell types.

    This class is intended to populate externally allocated arrays with content
    that is compatible with what VTK expects. This approach allows an improved
    separation of the OpenFOAM mesh description and the storage, and allows
    support of alternative storage containers (eg, std::vector, vtkDataArray).
    The ideal goal would be a zero-copy mechanism, but this does not work for
    several reasons:
    \par
    - OpenFOAM and VTK have different point ordering for prism
    - polyhedral decomposition
    - face-stream are required for VTK
    - VTK internal storage includes list size as part of the data
    - VTK includes storage may be a different base size (eg, long long)
      compared to the OpenFOAM label.

    \par Data Entries (slots)

    These are the storage entries normally associate with each output-type:
    \table
        legacy output (vtk DataFile Version 2.0)
        \c types    | vtk cell type (1-255)
        \c cells    | nLabels and unique vertex labels used by the cell, or
                    | [nLabels nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
    \endtable

    \table
        xml output
        \c types        | vtk cell type (1-255)
        \c connectivity | unique vertex labels used by the cell
        \c offsets      | end offset for each of \c connectivity
        \c faces        | face stream for polyhedral cells
                        | [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
        \c faceoffsets  | end offset for each of \c faces, with -1 for primitive cells
    \endtable

    \table
        internal1 storage (VTK-8 and earlier)
        \c types        | vtk cell type (1-255)
        \c connectivity | nLabels and unique vertex labels used by the cell
        \c location     | begin location for each of \c connectivity
        \c faces        | face stream for polyhedral cells
                        | [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
        \c facelocation | begin location for each of \c faces, with -1 for primitive cells
    \endtable

    \table
        internal2 storage (with VTK_CELL_ARRAY_V2)
        \c types        | vtk cell type (1-255)
        \c connectivity | unique vertex labels used by the cell
        \c offsets      | begin/end offsets for \c connectivity
        \c faces        | face stream for polyhedral cells
                        | [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
        \c facelocation | begin location for each of \c faces, with -1 for primitive cells
    \endtable

    The VTK storage concept for "connectivity" and "faces" somewhat resemble
    a CompactListList.

Note
    It is possible to specify a global point offset (via the globalIndex)
    so that the cell point labels will use global numbering.
    There is no support for point renumbering with merged mesh points,
    since it likely more efficient to use VTK point-blanking to mark duplicate
    points instead of merging points ourselves.

Note
    The VTK_CELL_ARRAY_V2 define (from vtkCellArray.h) indicates if the new
    (internal2) new format is being used.

SourceFiles
    foamVtuSizing.C
    foamVtuSizingI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtuSizing_H
#define Foam_vtk_vtuSizing_H

#include "label.H"
#include "labelList.H"
#include "foamVtkMeshMaps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellShape;
class polyMesh;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                    Class Foam::vtk::vtuSizing Declaration
\*---------------------------------------------------------------------------*/

class vtuSizing
{
public:

    // Public Data

    //- Types of content that the storage may represent
    enum contentType : char
    {
        LEGACY,     //!< Legacy VTK content
        XML,        //!< XML (VTU) content
        INTERNAL1,  //!< Internal vtkUnstructuredGrid content
        INTERNAL2   //!< Internal vtkUnstructuredGrid content, VTK_CELL_ARRAY_V2
    };

    //- The possible storage 'slots' that can be used
    enum slotType : char
    {
        CELLS,         //!< Cell connectivity (ALL)
        CELLS_OFFSETS, //!< Cell end-offsets (XML), locations (INTERNAL1)
                       //!< or begin/end offsets (INTERNAL2)
        FACES,         //!< Face-stream (XML, INTERNAL)
        FACES_OFFSETS  //!< Faces end-offsets (XML) or locations (INTERNAL1)
    };

    //- How the mesh cells have been selected or defined
    enum selectionModeType : char
    {
        FULL_MESH,
        SUBSET_MESH,
        SHAPE_MESH
    };


private:

    // Private Member Data

        //- Polyhedral decomposition requested
        bool decompose_;

        //- Manifold cells detected
        bool manifold_;

        //- How the mesh cells have been selected or defined
        selectionModeType selectionMode_;

        //- Number of cells in the mesh
        label nCells_;

        //- Number of points in the mesh
        label nPoints_;

        //- Number of vertex labels to represent the mesh
        label nVertLabels_;


    // Polyhedrals

        //- Number of polyhedral face labels for the mesh
        label nFaceLabels_;

        //- Number of polyhedral cells (informational)
        label nCellsPoly_;

        //- Number of vertex labels used by polyhedrals
        label nVertPoly_;


    // Decomposed polyhedrals

        //- Number of additional (decomposed) cells for the mesh
        label nAddCells_;

        //- Number of additional (decomposed) points for the mesh
        label nAddPoints_;

        //- Number of additional (decomposed) vertices for the mesh
        label nAddVerts_;


    // Private Member Functions

        //- set-size for cellMap and additionalIds
        void presizeMaps(foamVtkMeshMaps& maps) const;

        //- Verify list sizes
        static void checkSizes
        (
            const vtk::vtuSizing& sizing,

            const label cellTypes_size,
            const label vertLabels_size,
            const label vertOffset_size,
            const label faceLabels_size,
            const label faceOffset_size,

            const enum contentType output,

            const label cellMap_size,
            const label addPointsIds_size
        );

        //- Adjust cell/face offsets
        template<class LabelType>
        static void adjustOffsets
        (
            UList<LabelType>& vertOffset,
            UList<LabelType>& faceOffset,
            const enum contentType output,
            const bool hasFaceStream
        );

        //- Populate lists. For (legacy | xml | internal) VTK representations
        template<class LabelType>
        static void populateArrays
        (
            const polyMesh& mesh,
            const vtk::vtuSizing& sizing,
            UList<uint8_t>& cellTypes,
            UList<LabelType>& vertLabels,
            UList<LabelType>& vertOffset,
            UList<LabelType>& faceLabels,
            UList<LabelType>& faceOffset,
            const enum contentType output,
            labelUList& cellMap,
            labelUList& addPointsIds
        );

        //- Populate lists - Primitive mesh shapes only!!
        template<class LabelType>
        static void populateArrays
        (
            const UList<cellShape>& shapes,
            const vtk::vtuSizing& sizing,
            UList<uint8_t>& cellTypes,
            UList<LabelType>& vertLabels,
            UList<LabelType>& vertOffset,
            UList<LabelType>& faceLabels,    // unused
            UList<LabelType>& faceOffset,    // unused
            const enum contentType output,
            labelUList& cellMap,
            labelUList& addPointsIds  // unused
        );


protected:

    // Protected Member Functions

        //- Reset list for primitive shapes only (ADVANCED USAGE)
        void populateShapesLegacy
        (
            const UList<cellShape>& shapes,
            UList<uint8_t>& cellTypes,
            labelUList& connectivity,
            foamVtkMeshMaps& maps
        ) const;

        //- Reset list for primitive shapes only (ADVANCED USAGE)
        void populateShapesXml
        (
            const UList<cellShape>& shapes,
            UList<uint8_t>& cellTypes,
            labelUList& connectivity,
            labelUList& offsets,
            labelUList& faces,
            labelUList& facesOffsets,
            foamVtkMeshMaps& maps
        ) const;


public:

    // Constructors

        //- Default construct
        vtuSizing() noexcept;

        //- Construct sizing by analyzing the mesh.
        //  No polyhedral decomposition.
        explicit vtuSizing(const polyMesh& mesh);

        //- Construct sizing by analyzing the mesh.
        //  Optionally with polyhedral decomposition.
        vtuSizing(const polyMesh& mesh, const bool decompose);


    // Member Functions

    // Edit

        //- Reset sizing by analyzing the mesh.
        //  Optionally with polyhedral decomposition.
        void reset(const polyMesh& mesh, const bool decompose=false);

        //- Reset sizing by analyzing a subset of the mesh.
        //  A labelUList::null() selection corresponds to the entire mesh.
        //  Polyhedral decomposition is disabled for subsets.
        void reset
        (
            const polyMesh& mesh,
            const labelUList& subsetCellsIds,
            const bool decompose = false
        );

        //- Reset sizing using primitive shapes only (ADVANCED USAGE)
        //  Effectively removes any polyhedrals!
        void resetShapes(const UList<cellShape>& shapes);

        //- Reset all sizes to zero.
        void clear() noexcept;


    // Access

        //- Query the decompose flag (normally off)
        inline bool decompose() const noexcept;

        //- Manifold mesh cells detected? Globally consistent quantity.
        inline bool manifold() const noexcept;

        //- Query how the mesh cells have been selected or defined
        inline selectionModeType selectionMode() const noexcept;

        //- Number of cells for the mesh
        inline label nCells() const noexcept;

        //- Number of points for the mesh
        inline label nPoints() const noexcept;

        //- Number of vertex labels for the mesh
        inline label nVertLabels() const noexcept;

        //- Number of polyhedral face labels for the mesh
        inline label nFaceLabels() const noexcept;

        //- Number of polyhedral cells for the mesh
        inline label nCellsPoly() const noexcept;

        //- Number of vertex labels for polyhedral cells of the mesh
        inline label nVertPoly() const noexcept;

        //- Number of additional (decomposed) cells for the mesh
        inline label nAddCells() const noexcept;

        //- Number of additional (decomposed) points for the mesh
        inline label nAddPoints() const noexcept;

        //- Number of additional (decomposed) vertices for the mesh
        inline label nAddVerts() const noexcept;


        //- Number of field cells = nCells + nAddCells
        inline label nFieldCells() const noexcept;

        //- Number of field points = nPoints + nAddPoints
        inline label nFieldPoints() const noexcept;


    // Edit

        //- Alter number of mesh points (ADVANCED USAGE)
        void setNumPoints(label n) noexcept { nPoints_ = n; }

        //- Alter number of additional (cell-centre) points (ADVANCED USAGE)
        void setNumAddPoints(label n) noexcept { nAddPoints_ = n; }


    // Derived Sizes

        //- Return the required size for the storage slot
        label sizeOf
        (
            const enum contentType output,
            const enum slotType slot
        ) const;


        //- The calculated size for legacy storage
        inline label sizeLegacy() const;

        //- The calculated size for legacy storage of the specified slot
        inline label sizeLegacy(const enum slotType slot) const;

        //- The calculated size for xml storage of the specified slot
        inline label sizeXml(const enum slotType slot) const;

        //- The calculated size for vtk-internal storage of the specified slot
        inline label sizeInternal1(const enum slotType slot) const;

        //- The calculated size for vtk-internal storage of the specified slot
        inline label sizeInternal2(const enum slotType slot) const;


    // Routines for populating the output lists

        //- Populate lists for Legacy output
        void populateLegacy
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            labelUList& connectivity,
            foamVtkMeshMaps& maps
        ) const;

        //- Populate lists for XML output
        void populateXml
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            labelUList& connectivity,
            labelUList& offsets,
            labelUList& faces,
            labelUList& facesOffsets,
            foamVtkMeshMaps& maps
        ) const;


    // Internal types. The size of vtkIdType is unknown here

#undef  declarePopulateInternalMethod
#define declarePopulateInternalMethod(Type)                                   \
                                                                              \
        /*! \brief Populate lists for Internal VTK format */                  \
        void populateInternal                                                 \
        (                                                                     \
            const polyMesh& mesh,                                             \
            UList<uint8_t>& cellTypes,                                        \
            UList<Type>& connectivity,                                        \
            UList<Type>& offsets,                                             \
            UList<Type>& faces,                                               \
            UList<Type>& facesOffsets,                                        \
            foamVtkMeshMaps& maps,                                            \
            const enum contentType output                                     \
        ) const;                                                              \
                                                                              \
        /*! \brief Populate lists for Internal VTK format */                  \
        void populateInternal                                                 \
        (                                                                     \
            const polyMesh& mesh,                                             \
            UList<uint8_t>& cellTypes,                                        \
            UList<Type>& connectivity,                                        \
            UList<Type>& offsets,                                             \
            UList<Type>& faces,                                               \
            UList<Type>& facesOffsets,                                        \
            labelUList& cellMap,                                              \
            labelUList& addPointsIds,                                         \
            const enum contentType output                                     \
        ) const


        declarePopulateInternalMethod(int);
        declarePopulateInternalMethod(long);
        declarePopulateInternalMethod(long long);

        #undef declarePopulateInternalMethod


    // Routines for renumber vertices with a global point offset
    // Legacy and xml only, internal version less likely to be needed

        //- Copy vertex labels with a global point offset - legacy format
        static labelList copyVertLabelsLegacy
        (
            const labelUList& connectivity,
            const label globalPointOffset
        );

        //- Copy vertex labels with a global point offset - XML format
        static labelList copyVertLabelsXml
        (
            const labelUList& connectivity,
            const label globalPointOffset
        );

        //- Copy faces stream labels with a global point offset - XML format
        static labelList copyFaceLabelsXml
        (
            const labelUList& faceLabels,
            const label globalPointOffset
        );

        //- Copy face offsets with an offset from previous - XML format
        static labelList copyFaceOffsetsXml
        (
            const labelUList& faceOffsets,
            const label prevOffset
        );

        //- Renumber vertex labels by global point offset - legacy format
        static void renumberVertLabelsLegacy
        (
            labelUList& connectivity,
            const label globalPointOffset
        );

        //- Renumber vertex labels by global point offset - XML format
        static void renumberVertLabelsXml
        (
            labelUList& connectivity,
            const label globalPointOffset
        );

        //- Renumber faces stream labels by global point offset - XML format
        static void renumberFaceLabelsXml
        (
            labelUList& faceLabels,
            const label globalPointOffset
        );

        //- Renumber face offsets with an offset from previous - XML format
        static void renumberFaceOffsetsXml
        (
            labelUList& faceOffsets,
            const label prevOffset
        );


    // Write

        //- Report some information
        void info(Ostream& os) const;


    // Member Operators

        //- Test equality
        bool operator==(const vtuSizing& rhs) const;

        //- Test inequality
        bool operator!=(const vtuSizing& rhs) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtuSizingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
