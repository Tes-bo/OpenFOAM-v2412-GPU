/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2021 OpenCFD Ltd.
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
    Foam::vtk::vtuCells

Description
    A deep-copy description of an OpenFOAM volume mesh in data structures
    suitable for VTK UnstructuredGrid, including the possibility of
    decomposing polyhedral cells into primitive cell types.

    Knowledge of the vtkUnstructuredGrid and the corresponding \c .vtu
    xml file-format aids in understanding this class.
    The class can be used for the VTK xml format, legacy format, as well as a
    VTK internal representation. The internal representation is somewhat
    related to the xml format, but not entirely.

SeeAlso
    Foam::vtk::vtuSizing

SourceFiles
    foamVtuCells.C
    foamVtuCellsI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtuCells_H
#define Foam_vtk_vtuCells_H

#include "foamVtkCore.H"
#include "foamVtkMeshMaps.H"
#include "foamVtuSizing.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace vtk
{
// Forward Declarations
class outputOptions;

/*---------------------------------------------------------------------------*\
                     Class Foam::vtk::vtuCells Declaration
\*---------------------------------------------------------------------------*/

class vtuCells
:
    public vtuSizing
{
    // Private Data

    // Requested output types

        //- Output content type
        contentType output_;

        //- Bookkeeping for polyhedral cell decomposition
        bool decomposeRequest_;


    // Storage of output

        //- Cell types (including added cells) in vtk numbering
        //  Range is 1-255
        List<uint8_t> cellTypes_;

        //- Vertices per cell (including added cells) in vtk ordering
        List<label> vertLabels_;

        //- Connectivity (vertices) offset for the end of each cell
        List<label> vertOffset_;

        //- Face lists per polyhedral cell
        List<label> faceLabels_;

        //- Face label offsets
        List<label> faceOffset_;

        //- Bookkeeping for polyhedral cell decomposition
        foamVtkMeshMaps maps_;


    // Private Member Functions

        //- Resize storage
        void resize_all();

        //- Create the geometry entries
        void populateOutput(const polyMesh& mesh);

        //- Create the geometry entries from shapes
        void populateOutput(const UList<cellShape>& shapes);

        //- No copy construct
        vtuCells(const vtuCells&) = delete;

        //- No copy assignment
        void operator=(const vtuCells&) = delete;


public:

    // Constructors

        //- Default construct (XML format, no polyhedral decomposition)
        explicit vtuCells
        (
            const enum contentType output = contentType::XML,
            const bool decompose = false
        );

        //- Construct from components, create output information immediately
        explicit vtuCells
        (
            const polyMesh& mesh,
            const enum contentType output = contentType::XML,
            const bool decompose = false
        );

        //- Construct from components.
        //  Optionally with polyhedral decomposition.
        explicit vtuCells
        (
            const vtk::outputOptions opts,
            const bool decompose = false
        );

        //- Construct from components, create output information immediately
        vtuCells
        (
            const polyMesh& mesh,
            const vtk::outputOptions opts,
            const bool decompose = false
        );


    // Member Functions

    // Access

        //- The output content type
        inline enum contentType content() const noexcept;

        //- Query the polyhedral decompose requested flag.
        inline bool decomposeRequested() const noexcept;

        //- True if no cellTypes are populated.
        inline bool empty() const noexcept;

        //- The size of populated cellTypes.
        inline label size() const noexcept;


    // Edit

        //- Reset all sizes to zero.
        void clear();

        //- Create the geometry using the previously requested output and
        //- decomposition types.
        void reset(const polyMesh& mesh);

        //- Create the geometry for a mesh subset,
        //- using previously requested output and decomposition types.
        void reset
        (
            const polyMesh& mesh,
            const labelUList& subsetCellsIds
        );

        //- Respecify requested output and decomposition type prior to
        //- creating the geometry
        void reset
        (
            const polyMesh& mesh,
            const enum contentType output,
            const bool decompose
        );

        //- Reset sizing using primitive shapes only (ADVANCED USAGE)
        //  Effectively removes any polyhedrals!
        void resetShapes(const UList<cellShape>& shapes);

        //- Renumber cell ids to account for subset meshes
        void renumberCells(const labelUList& mapping);

        //- Renumber point ids to account for subset meshes
        void renumberPoints(const labelUList& mapping);

        //- Define which additional cell-centres are to be used (ADVANCED!)
        void addPointCellLabels(const labelUList& cellIds);


    // Storage Access

        //- Values for "types" (XML) and "CELL_TYPES" (legacy)
        inline const List<uint8_t>& cellTypes() const noexcept;

        //- Values for "connectivity" (XML) or "CELLS" (legacy)
        inline const labelList& vertLabels() const noexcept;

        //- Values for "offsets" (XML only)
        inline const labelList& vertOffsets() const noexcept;

        //- Values for "faces" (XML only)
        inline const labelList& faceLabels() const noexcept;

        //- Values for "faceoffset" (XML only)
        inline const labelList& faceOffsets() const noexcept;

        //- Additional point addressing (from added point to original cell)
        inline const labelList& addPointCellLabels() const noexcept;

        //- Original cell ids for all cells (regular and decomposed).
        inline const labelList& cellMap() const noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtuCellsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
