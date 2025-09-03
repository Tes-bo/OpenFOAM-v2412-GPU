/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "foamVtuSizing.H"
#include "foamVtkCore.H"
#include "polyMesh.H"
#include "cellShape.H"
#include "manifoldCellsMeshObject.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class LabelType>
void Foam::vtk::vtuSizing::adjustOffsets
(
    UList<LabelType>& vertOffset,
    UList<LabelType>& faceOffset,
    const enum contentType output,
    const bool hasFaceStream
)
{
    // ===========================================
    // Adjust vertOffset for all cells
    // A second pass is needed for several reasons.
    // - Additional (decomposed) cells are placed out of sequence
    // - INTERNAL1 connectivity has size prefixed
    //
    // Cell offsets:
    // - XML format expects end-offsets,
    // - INTERNAL1 expects begin-offsets
    // - INTERNAL2 expects begin/end-offsets

    switch (output)
    {
        case contentType::LEGACY: // Nothing to do
            break;

        case contentType::XML:
        {
            // Transform cell sizes (vertOffset) into begin offsets

            // vertOffset[0] already contains its size, leave untouched
            for (label i = 1; i < vertOffset.size(); ++i)
            {
                vertOffset[i] += vertOffset[i-1];
            }

            // The end face offsets, leaving -1 untouched
            if (hasFaceStream)
            {
                LabelType prev(0);

                for (LabelType& off : faceOffset)
                {
                    const LabelType sz(off);
                    if (sz > 0)
                    {
                        prev += sz;
                        off = prev;
                    }
                }
            }
            break;
        }

        case contentType::INTERNAL1:
        {
            // Transform cell sizes (vertOffset) into begin offsets
            {
                LabelType beg(0);

                for (LabelType& off : vertOffset)
                {
                    const LabelType sz(off);
                    off = beg;
                    beg += 1 + sz;  // Additional 1 to skip embedded prefix
                }
            }

            // The begin face offsets, leaving -1 untouched
            if (hasFaceStream)
            {
                LabelType beg(0);

                for (LabelType& off : faceOffset)
                {
                    const LabelType sz(off);
                    if (sz > 0)
                    {
                        off = beg;
                        beg += sz;
                    }
                }
            }
            break;
        }

        case contentType::INTERNAL2:
        {
            // Transform cell sizes (vertOffset) into begin/end offsets
            // input    [n1, n2, n3, ..., 0]
            // becomes  [0, n1, n1+n2, n1+n2+n3, ..., nTotal]

            // The last entry of vertOffset was initialized as zero and
            // never revisited, so the following loop is OK
            {
                LabelType total(0);

                for (LabelType& off : vertOffset)
                {
                    const LabelType sz(off);
                    off = total;
                    total += sz;
                }
            }

            // The begin face offsets, leaving -1 untouched
            if (hasFaceStream)
            {
                LabelType beg(0);

                for (LabelType& off : faceOffset)
                {
                    const LabelType sz(off);
                    if (sz > 0)
                    {
                        off = beg;
                        beg += sz;
                    }
                }
            }
            break;
        }
    }
}


template<class LabelType>
void Foam::vtk::vtuSizing::populateArrays
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
)
{
    if (sizing.selectionMode() == selectionModeType::SHAPE_MESH)
    {
        FatalErrorInFunction
            << "Programming error ... attempting to populate a VTU mesh"
            << " but it was originally sized using independent cell shapes"
            << exit(FatalError);
    }

    // Verify storage sizes
    checkSizes
    (
        sizing,

        cellTypes.size(),
        vertLabels.size(), vertOffset.size(),
        faceLabels.size(), faceOffset.size(),

        output,

        cellMap.size(),
        addPointsIds.size()
    );

    // Characteristics

    // Are vertLabels prefixed with the size?
    // Also use as the size of the prefixed information
    const int prefix =
    (
        output == contentType::LEGACY
     || output == contentType::INTERNAL1
    ) ? 1 : 0;


    // Initialization

    faceOffset = -1;

    // For INTERNAL2, the vertOffset is (nFieldCells+1), which means that
    // the last entry is never visited. Set as zero now.

    if (vertOffset.size())
    {
        vertOffset.first() = 0;
        vertOffset.last() = 0;
    }


    const cellModel& tet      = cellModel::ref(cellModel::TET);
    const cellModel& pyr      = cellModel::ref(cellModel::PYR);
    const cellModel& prism    = cellModel::ref(cellModel::PRISM);
    const cellModel& hex      = cellModel::ref(cellModel::HEX);
    const cellModel& wedge    = cellModel::ref(cellModel::WEDGE);
    const cellModel& tetWedge = cellModel::ref(cellModel::TETWEDGE);

    const cellShapeList& shapes = mesh.cellShapes();
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();

    // The face owner is needed to determine the face orientation
    const labelList& owner = mesh.faceOwner();

    // Unique vertex labels per polyhedral
    labelHashSet hashUniqId(512);

    // Index into vertLabels, faceLabels for normal cells
    label nVertLabels = 0;
    label nFaceLabels = 0;

    // Index into vertLabels for decomposed polys
    label nVertDecomp = sizing.nVertLabels() + prefix*sizing.nCells();

    // Placement of additional decomposed cells
    label nCellDecomp = mesh.nCells();

    // Placement of additional point labels
    label nPointDecomp = mesh.nPoints();

    // Non-decomposed polyhedral are represented as a face-stream.
    // For legacy format, this stream replaces the normal connectivity
    // information. Use references to alias where the face output should land.

    UList<LabelType>& faceOutput =
    (
        output == contentType::LEGACY
      ? vertLabels
      : faceLabels
    );

    label& faceIndexer =
    (
        output == contentType::LEGACY
      ? nVertLabels
      : nFaceLabels
    );

    // ===========================================
    // STAGE 2: Rewrite in VTK form
    // During this stage, the vertOffset contains the *size* associated with
    // the per-cell vertLabels entries, and the faceOffset contains the *size*
    // associated with the per-cell faceLabels.


    // Special treatment for mesh subsets
    // Here the cellMap is the list of input cells!

    const bool isSubsetMesh
    (
        sizing.selectionMode() == selectionModeType::SUBSET_MESH
    );

    const label nInputCells =
    (
        isSubsetMesh
      ? cellMap.size()
      : shapes.size()
    );


    for
    (
        label inputi = 0, cellIndex = 0; // cellIndex: the ouput location
        inputi < nInputCells;
        ++inputi, ++cellIndex
    )
    {
        const label celli(isSubsetMesh ? cellMap[inputi] : inputi);

        const cellShape& shape = shapes[celli];
        const cellModel& model = shape.model();

        if (!isSubsetMesh)
        {
            cellMap[cellIndex] = celli;
        }

        if (model == tet)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_TETRA;
            constexpr label nShapePoints = 4; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == pyr)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_PYRAMID;
            constexpr label nShapePoints = 5; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == hex)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_HEXAHEDRON;
            constexpr label nShapePoints = 8; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == prism)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_WEDGE;
            constexpr label nShapePoints = 6; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            // VTK_WEDGE triangles point outwards (swap 1<->2, 4<->5)
            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[4];
        }
        else if (model == tetWedge && sizing.decompose())
        {
            // Treat as squeezed prism
            cellTypes[cellIndex] = vtk::cellType::VTK_WEDGE;
            constexpr label nShapePoints = 6;

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[4];
            vertLabels[nVertLabels++] = shape[3];
        }
        else if (model == wedge && sizing.decompose())
        {
            // Treat as squeezed hex
            cellTypes[cellIndex] = vtk::cellType::VTK_HEXAHEDRON;
            constexpr label nShapePoints = 8;

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[4];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[6];
        }
        else if (sizing.decompose())
        {
            // Polyhedral cell - decompose into tet/pyr.

            // Ensure we have the correct orientation for the base of the
            // primitive cell shape.
            // If the cell is face owner, the orientation needs to be flipped
            // to avoid defining negative cells.
            // VTK may not care, but we'll do it anyhow for safety.

            // Mapping from additional point to cell, and the new vertex from
            // the cell-centre
            const label newVertexLabel = nPointDecomp;

            addPointsIds[nPointDecomp++] = celli;

            // Whether to insert cell in place of original or not.
            bool firstCell = true;

            const labelList& cFaces = meshCells[celli];

            for (const label facei : cFaces)
            {
                const face& f = meshFaces[facei];
                const bool isOwner = (owner[facei] == celli);

                // Count triangles/quads in decomposition
                label nTria = 0, nQuad = 0;
                f.nTrianglesQuads(mesh.points(), nTria, nQuad);

                // Do actual decomposition
                faceList faces3(nTria);
                faceList faces4(nQuad);
                nTria = 0, nQuad = 0;
                f.trianglesQuads(mesh.points(), nTria, nQuad, faces3, faces4);

                for (const face& quad : faces4)
                {
                    // Quad becomes a pyramid

                    constexpr label nShapePoints = 5;  // pyr (5 vertices)

                    label celLoc, vrtLoc;
                    if (firstCell)
                    {
                        firstCell = false;
                        celLoc = cellIndex;
                        vrtLoc = nVertLabels;
                        nVertLabels += prefix + nShapePoints;
                    }
                    else
                    {
                        celLoc = nCellDecomp++;
                        vrtLoc = nVertDecomp;
                        nVertDecomp += prefix + nShapePoints;
                    }
                    cellMap[celLoc] = celli;

                    cellTypes[celLoc] = vtk::cellType::VTK_PYRAMID;
                    if (vertOffset.size())
                    {
                        vertOffset[celLoc] = nShapePoints;
                    }
                    if (prefix)
                    {
                        vertLabels[vrtLoc++] = nShapePoints;
                    }

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels[vrtLoc++] = quad[0];
                        vertLabels[vrtLoc++] = quad[3];
                        vertLabels[vrtLoc++] = quad[2];
                        vertLabels[vrtLoc++] = quad[1];
                    }
                    else
                    {
                        vertLabels[vrtLoc++] = quad[0];
                        vertLabels[vrtLoc++] = quad[1];
                        vertLabels[vrtLoc++] = quad[2];
                        vertLabels[vrtLoc++] = quad[3];
                    }

                    // The apex
                    vertLabels[vrtLoc++] = newVertexLabel;
                }

                for (const face& tria : faces3)
                {
                    // Triangle becomes a tetrahedral

                    constexpr label nShapePoints = 4;  // tet (4 vertices)

                    label celLoc, vrtLoc;
                    if (firstCell)
                    {
                        firstCell = false;
                        celLoc = cellIndex;
                        vrtLoc = nVertLabels;
                        nVertLabels += prefix + nShapePoints;
                    }
                    else
                    {
                        celLoc = nCellDecomp++;
                        vrtLoc = nVertDecomp;
                        nVertDecomp += prefix + nShapePoints;
                    }
                    cellMap[celLoc] = celli;

                    cellTypes[celLoc] = vtk::cellType::VTK_TETRA;
                    if (vertOffset.size())
                    {
                        vertOffset[celLoc] = nShapePoints;
                    }
                    if (prefix)
                    {
                        vertLabels[vrtLoc++] = nShapePoints;
                    }

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels[vrtLoc++] = tria[0];
                        vertLabels[vrtLoc++] = tria[2];
                        vertLabels[vrtLoc++] = tria[1];
                    }
                    else
                    {
                        vertLabels[vrtLoc++] = tria[0];
                        vertLabels[vrtLoc++] = tria[1];
                        vertLabels[vrtLoc++] = tria[2];
                    }

                    // The apex
                    vertLabels[vrtLoc++] = newVertexLabel;
                }
            }
        }
        else
        {
            // Polyhedral cell - not decomposed
            hashUniqId.clear();  // unique node ids used (XML, INTERNAL)

            // face-stream
            //   [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
            cellTypes[cellIndex] = vtk::cellType::VTK_POLYHEDRON;

            const labelList& cFaces = meshCells[celli];

            const label startLabel = faceIndexer;

            if (output == contentType::LEGACY)
            {
                faceOutput[startLabel] = 0; // placeholder for total size
                ++faceIndexer;
            }

            faceOutput[faceIndexer++] = cFaces.size();

            for (const label facei : cFaces)
            {
                const face& f = mesh.faces()[facei];
                const bool isOwner = (owner[facei] == celli);
                const label nFacePoints = f.size();

                hashUniqId.insert(f);

                // The number of labels for this face
                faceOutput[faceIndexer++] = nFacePoints;

                faceOutput[faceIndexer++] = f[0];
                if (isOwner)
                {
                    for (label fp = 1; fp < nFacePoints; ++fp)
                    {
                        faceOutput[faceIndexer++] = f[fp];
                    }
                }
                else
                {
                    for (label fp = nFacePoints - 1; fp > 0; --fp)
                    {
                        faceOutput[faceIndexer++] = f[fp];
                    }
                }
            }

            if (output == contentType::LEGACY)
            {
                // Update size for legacy face stream
                // (subtract 1 to avoid counting the storage location)
                faceOutput[startLabel] = (faceIndexer - 1 - startLabel);
            }
            else
            {
                // Size for face stream
                faceOffset[cellIndex] = (faceIndexer - startLabel);

                vertOffset[cellIndex] = hashUniqId.size();
                if (prefix)
                {
                    vertLabels[nVertLabels++] = hashUniqId.size();
                }

                for (const label pointi : hashUniqId.sortedToc())
                {
                    vertLabels[nVertLabels++] = pointi;
                }
            }
        }
    }

    // ===========================================
    // STAGE 3: Adjust vertOffset for all cells
    // A second pass is needed for several reasons.
    // - Additional (decomposed) cells are placed out of sequence
    // - INTERNAL1 connectivity has size prefixed
    //
    // Cell offsets:
    // - XML format expects end-offsets,
    // - INTERNAL1 expects begin-offsets
    // - INTERNAL2 expects begin/end-offsets

    adjustOffsets<LabelType>
    (
        vertOffset,
        faceOffset,
        output,
        sizing.nFaceLabels()  // hasFaceStream
    );
}



// Synchronize changes here with the following:
// - vtuSizing::resetShapes
// - vtuSizing::populateArrays

template<class LabelType>
void Foam::vtk::vtuSizing::populateArrays
(
    const UList<cellShape>& shapes,
    const vtk::vtuSizing& sizing,

    UList<uint8_t>& cellTypes,
    UList<LabelType>& vertLabels,
    UList<LabelType>& vertOffset,
    UList<LabelType>& faceLabels,
    UList<LabelType>& faceOffset,
    const enum contentType output,
    labelUList& cellMap,
    labelUList& addPointsIds
)
{
    if (sizing.selectionMode() != selectionModeType::SHAPE_MESH)
    {
        FatalErrorInFunction
            << "Programming error ... attempting to populate a VTU mesh"
            << " from cell shapes, but sizing originated from a different"
            << " representation" << nl
            << exit(FatalError);
    }

    // Verify storage sizes
    checkSizes
    (
        sizing,

        cellTypes.size(),
        vertLabels.size(), vertOffset.size(),
        faceLabels.size(), faceOffset.size(),

        output,

        cellMap.size(),
        addPointsIds.size()
    );

    // Characteristics

    // Are vertLabels prefixed with the size?
    // Also use as the size of the prefixed information
    const int prefix =
    (
        output == contentType::LEGACY
     || output == contentType::INTERNAL1
    ) ? 1 : 0;


    // Initialization

    faceOffset = -1;

    // For INTERNAL2, the vertOffset is (nFieldCells+1), which means that
    // the last entry is never visited. Set as zero now.

    if (vertOffset.size())
    {
        vertOffset.first() = 0;
        vertOffset.last() = 0;
    }


    const cellModel& tet      = cellModel::ref(cellModel::TET);
    const cellModel& pyr      = cellModel::ref(cellModel::PYR);
    const cellModel& prism    = cellModel::ref(cellModel::PRISM);
    const cellModel& hex      = cellModel::ref(cellModel::HEX);

    // Index into vertLabels for normal cells
    label nVertLabels = 0;

    // ===========================================
    // STAGE 2: Rewrite in VTK form
    // During this stage, the vertOffset contains the *size* associated with
    // the per-cell vertLabels entries, and the faceOffset contains the *size*
    // associated with the per-cell faceLabels.

    const label nInputCells = shapes.size();

    // label nIgnored = 0;

    for
    (
        label inputi = 0, cellIndex = 0; // cellIndex: the ouput location
        inputi < nInputCells;
        ++inputi, ++cellIndex
    )
    {
        const cellShape& shape = shapes[inputi];
        const cellModel& model = shape.model();

        if (model == tet)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_TETRA;
            constexpr label nShapePoints = 4; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == pyr)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_PYRAMID;
            constexpr label nShapePoints = 5; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == hex)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_HEXAHEDRON;
            constexpr label nShapePoints = 8; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == prism)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_WEDGE;
            constexpr label nShapePoints = 6; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            // VTK_WEDGE triangles point outwards (swap 1<->2, 4<->5)
            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[4];
        }
        else
        {
            // Silent here.
            // - already complained (and skipped) during initial sizing
            --cellIndex;
            // ++nIgnored;
        }
    }

    // May have been done by caller,
    // but for additional safety set an identity mapping
    Foam::identity(cellMap);

    // ===========================================
    // Adjust vertOffset for all cells
    // A second pass is needed for several reasons.
    // - Additional (decomposed) cells are placed out of sequence
    // - INTERNAL1 connectivity has size prefixed
    //
    // Cell offsets:
    // - XML format expects end-offsets,
    // - INTERNAL1 expects begin-offsets
    // - INTERNAL2 expects begin/end-offsets

    adjustOffsets<LabelType>
    (
        vertOffset,
        faceOffset,
        output,
        sizing.nFaceLabels()  // hasFaceStream
    );
}


//unused template<class LabelType, class LabelType2>
//unused void Foam::vtk::vtuSizing::renumberVertLabelsInternalImpl
//unused (
//unused     UList<uint8_t>& cellTypes,
//unused     UList<LabelType>& vertLabels,
//unused     const LabelType2 globalPointOffset
//unused )
//unused {
//unused     // INTERNAL vertLabels = "connectivity" contain
//unused     // [nLabels, vertex labels...]
//unused
//unused     auto iter = vertLabels.begin();
//unused     const auto last = vertLabels.end();
//unused
//unused     while (iter < last)
//unused     {
//unused         LabelType nLabels = *iter;
//unused         ++iter;
//unused
//unused         while (nLabels--)
//unused         {
//unused             *iter += globalPointOffset;
//unused             ++iter;
//unused         }
//unused     }
//unused }


// ************************************************************************* //
