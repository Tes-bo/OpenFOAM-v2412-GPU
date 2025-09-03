/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

// Only used in this file
#include "foamVtuSizingImpl.C"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::vtuSizing::presizeMaps(foamVtkMeshMaps& maps) const
{
    maps.cellMap().resize(this->nFieldCells());
    maps.additionalIds().resize(this->nAddPoints());
}


void Foam::vtk::vtuSizing::checkSizes
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
)
{
    label nErrors = 0;

    #undef  CHECK_SIZING
    #define CHECK_SIZING(what, sizeInput, sizeExpected)        \
    if (sizeInput != sizeExpected)                             \
    {                                                          \
        if (!nErrors++)                                        \
        {                                                      \
            FatalErrorInFunction << "VTK sizing error" << nl;  \
        }                                                      \
        FatalError                                             \
            << "    " << what << " size=" << sizeInput         \
            << " expected " << sizeExpected << nl;             \
    }


    CHECK_SIZING("cellTypes", cellTypes_size, sizing.nFieldCells());
    CHECK_SIZING("cellMap", cellMap_size, sizing.nFieldCells());
    CHECK_SIZING("addPointsIds", addPointsIds_size, sizing.nAddPoints());

    switch (output)
    {
        case contentType::LEGACY:
        {
            CHECK_SIZING("legacy", vertLabels_size, sizing.sizeLegacy());
            break;
        }

        case contentType::XML:
        {
            // XML uses connectivity/offset pair.
            CHECK_SIZING
            (
                "connectivity",
                vertLabels_size,
                sizing.sizeXml(slotType::CELLS)
            );
            CHECK_SIZING
            (
                "offsets",
                vertOffset_size,
                sizing.sizeXml(slotType::CELLS_OFFSETS)
            );
            if (sizing.nFaceLabels())
            {
                CHECK_SIZING
                (
                    "faces",
                    faceLabels_size,
                    sizing.sizeXml(slotType::FACES)
                );

                CHECK_SIZING
                (
                    "faceOffsets",
                    faceOffset_size,
                    sizing.sizeXml(slotType::FACES_OFFSETS)
                );
            }
            break;
        }

        case contentType::INTERNAL1:
        {
            // VTK-internal1 connectivity/offset pair.
            CHECK_SIZING
            (
                "connectivity",
                vertLabels_size,
                sizing.sizeInternal1(slotType::CELLS)
            );
            CHECK_SIZING
            (
                "offsets",
                vertOffset_size,
                sizing.sizeInternal1(slotType::CELLS_OFFSETS)
            );
            if (sizing.nFaceLabels())
            {
                CHECK_SIZING
                (
                    "faces",
                    faceLabels_size,
                    sizing.sizeInternal1(slotType::FACES)
                );
                CHECK_SIZING
                (
                    "faceOffsets",
                    faceOffset_size,
                    sizing.sizeInternal1(slotType::FACES_OFFSETS)
                );
            }
            break;
        }

        case contentType::INTERNAL2:
        {
            // VTK-internal2 connectivity/offset pair.
            CHECK_SIZING
            (
                "connectivity",
                vertLabels_size,
                sizing.sizeInternal2(slotType::CELLS)
            );
            CHECK_SIZING
            (
                "offsets",
                vertOffset_size,
                sizing.sizeInternal2(slotType::CELLS_OFFSETS)
            );
            if (sizing.nFaceLabels())
            {
                CHECK_SIZING
                (
                    "faces",
                    faceLabels_size,
                    sizing.sizeInternal2(slotType::FACES)
                );
                CHECK_SIZING
                (
                    "faceOffsets",
                    faceOffset_size,
                    sizing.sizeInternal2(slotType::FACES_OFFSETS)
                );
            }
            break;
        }
    }

    if (nErrors)
    {
        FatalError
            << nl
            << "Total of " << nErrors << " sizing errors encountered!"
            << exit(FatalError);
    }

    #undef CHECK_SIZING
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::vtuSizing::vtuSizing() noexcept
{
    clear();
}


Foam::vtk::vtuSizing::vtuSizing
(
    const polyMesh& mesh,
    const bool decompose
)
{
    clear();
    reset(mesh, decompose);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::vtuSizing::clear() noexcept
{
    decompose_   = false;
    manifold_    = false;
    selectionMode_ = FULL_MESH;
    nCells_      = 0;
    nPoints_     = 0;
    nVertLabels_ = 0;

    nFaceLabels_ = 0;
    nCellsPoly_  = 0;
    nVertPoly_   = 0;

    nAddCells_   = 0;
    nAddPoints_  = 0;
    nAddVerts_   = 0;
}


void Foam::vtk::vtuSizing::reset
(
    const polyMesh& mesh,
    const bool decompose
)
{
    reset(mesh, labelUList::null(), decompose);
}


void Foam::vtk::vtuSizing::reset
(
    const polyMesh& mesh,
    const labelUList& subsetCellsIds,
    const bool decompose
)
{
    // References to cell shape models
    const cellModel& tet   = cellModel::ref(cellModel::TET);
    const cellModel& pyr   = cellModel::ref(cellModel::PYR);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& hex   = cellModel::ref(cellModel::HEX);
    const cellModel& wedge    = cellModel::ref(cellModel::WEDGE);
    const cellModel& tetWedge = cellModel::ref(cellModel::TETWEDGE);

    const cellShapeList& shapes = mesh.cellShapes();
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();

    // Unique vertex labels per polyhedral
    labelHashSet hashUniqId(2*256);


    // Special treatment for mesh subsets.
    const bool isSubsetMesh
    (
        notNull(subsetCellsIds)
    );

    if (isSubsetMesh)
    {
        decompose_  = false;  // Disallow decomposition for subset mode
        selectionMode_ = selectionModeType::SUBSET_MESH;
    }
    else
    {
        decompose_  = decompose;  // Disallow decomposition
        selectionMode_ = selectionModeType::FULL_MESH;
    }

    // Manifold cells detected?
    manifold_ = manifoldCellsMeshObject::New(mesh).manifold();

    const label nInputCells =
    (
        isSubsetMesh
      ? subsetCellsIds.size()
      : shapes.size()
    );

    nCells_    = nInputCells;
    nPoints_   = mesh.nPoints();
    nAddCells_ = 0;
    nAddVerts_ = 0;

    nCellsPoly_  = nCells_;
    nVertLabels_ = 0;
    nFaceLabels_ = 0;
    nVertPoly_   = 0;

    for (label inputi = 0; inputi < nInputCells; ++inputi)
    {
        const label celli(isSubsetMesh ? subsetCellsIds[inputi] : inputi);

        const cellShape& shape = shapes[celli];
        const cellModel& model = shape.model();

        if
        (
            model == tet
         || model == pyr
         || model == prism
         || model == hex
        )
        {
            // Normal primitive - not a poly
            --nCellsPoly_;
            nVertLabels_ += shape.size();
        }
        else if (model == tetWedge && decompose_)
        {
            nVertLabels_ += 6;  // Treat as squeezed prism (VTK_WEDGE)
        }
        else if (model == wedge && decompose_)
        {
            nVertLabels_ += 8;  // Treat as squeezed hex
        }
        else if (decompose_)
        {
            // Polyhedral: Decompose into tets + pyramids.
            ++nAddPoints_;

            // Count vertices into first decomposed cell
            bool first = true;

            const cell& cFaces = meshCells[celli];
            for (const label facei : cFaces)
            {
                const face& f = meshFaces[facei];

                // Face decomposed into triangles and quads
                // Tri -> Tet, Quad -> Pyr
                label nTria = 0, nQuad = 0;
                f.nTrianglesQuads(mesh.points(), nTria, nQuad);

                nAddCells_ += nTria + nQuad;
                nAddVerts_ += (nTria * 4) + (nQuad * 5);

                if (first)
                {
                    first = false;
                    --nAddCells_;

                    const label nvrt = (nQuad ? 5 : 4);
                    nAddVerts_   -= nvrt;
                    nVertLabels_ += nvrt;
                }
            }
        }
        else
        {
            // Polyhedral: Not decomposed

            const labelList& cFaces = meshCells[celli];

            // Unique node ids used (XML/INTERNAL, not needed for LEGACY)
            hashUniqId.clear();

            // Face stream sizing:
            // number of faces, size of each face, vertices per face
            // [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]

            for (const label facei : cFaces)
            {
                const face& f = meshFaces[facei];
                nFaceLabels_ += f.size();

                hashUniqId.insert(f);
            }

            // Legacy format only uses the face-stream.
            // - track what *NOT* to use for legacy
            nVertLabels_ += hashUniqId.size();
            nVertPoly_   += hashUniqId.size();

            nFaceLabels_ += 1 + cFaces.size();
        }
    }

    // Requested and actually required
    decompose_ = (decompose_ && nCellsPoly_);
}


// Synchronize changes here with the following:
// - vtuSizing::resetShapes
// - vtuSizing::populateArrays
//
void Foam::vtk::vtuSizing::resetShapes
(
    const UList<cellShape>& shapes
)
{
    const cellModel& tet      = cellModel::ref(cellModel::TET);
    const cellModel& pyr      = cellModel::ref(cellModel::PYR);
    const cellModel& prism    = cellModel::ref(cellModel::PRISM);
    const cellModel& hex      = cellModel::ref(cellModel::HEX);

    decompose_ = false;  // Disallow decomposition
    manifold_ = false;   // Assume no manifold cells possible

    selectionMode_ = SHAPE_MESH;

    const label nInputCells = shapes.size();

    nCells_    = nInputCells;
    nPoints_   = 0;
    nAddCells_ = 0;
    nAddVerts_ = 0;

    nCellsPoly_  = 0;
    nVertLabels_ = 0;
    nFaceLabels_ = 0;
    nVertPoly_   = 0;

    label nIgnored = 0;

    for (label inputi = 0; inputi < nInputCells; ++inputi)
    {
        const cellShape& shape = shapes[inputi];
        const cellModel& model = shape.model();

        if
        (
            model == tet
         || model == pyr
         || model == prism
         || model == hex
        )
        {
            nVertLabels_ += shape.size();

            // Guess for number of addressed points
            nPoints_ = max(nPoints_, max(shape));
        }
        else
        {
            --nCells_;
            ++nIgnored;
        }
    }

    if (nIgnored)
    {
        FatalErrorInFunction
            << "Encountered " << nIgnored << " unsupported cell shapes"
            << " ... this is likely not good" << nl
            << exit(FatalError);
    }

    if (nCells_)
    {
        ++nPoints_;
    }
}


Foam::label Foam::vtk::vtuSizing::sizeOf
(
    const enum contentType output,
    const enum slotType slot
) const
{
    switch (output)
    {
        case contentType::LEGACY:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // legacy uses connectivity for primitives, but directly
                    // stores face streams into connectivity as well.
                    // size-prefix per cell
                    return
                    (
                        nVertLabels() + nAddVerts() - nVertPoly() // primitives
                      + nFaceLabels()     // face-stream (poly)
                      + nFieldCells()     // nFieldCells (size prefix)
                    );
                    break;

                default:
                    break;
            }
            break;
        }

        case contentType::XML:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    return (nVertLabels() + nAddVerts());
                    break;

                case slotType::CELLS_OFFSETS:
                    return nFieldCells();
                    break;

                case slotType::FACES:
                    return nFaceLabels();
                    break;

                case slotType::FACES_OFFSETS:
                    return nFaceLabels() ? nFieldCells() : 0;
                    break;
            }
            break;
        }

        case contentType::INTERNAL1:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // size-prefix per cell
                    return (nVertLabels() + nAddVerts() + nFieldCells());
                    break;

                case slotType::CELLS_OFFSETS:
                    return nFieldCells();
                    break;

                case slotType::FACES:
                    return nFaceLabels();
                    break;

                case slotType::FACES_OFFSETS:
                    return nFaceLabels() ? nFieldCells() : 0;
                    break;
            }
            break;
        }

        case contentType::INTERNAL2:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    return (nVertLabels() + nAddVerts());
                    break;

                case slotType::CELLS_OFFSETS:
                    return (nFieldCells() + 1);
                    break;

                case slotType::FACES:
                    return nFaceLabels();
                    break;

                case slotType::FACES_OFFSETS:
                    return nFaceLabels() ? nFieldCells() : 0;
                    break;
            }
            break;
        }
    }

    return 0;
}


// * * * * * * * * * * * * * *  Populate Lists * * * * * * * * * * * * * * * //

void Foam::vtk::vtuSizing::populateLegacy
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    labelUList& vertLabels,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        vertLabels,
        unused, // offsets
        unused, // faces
        unused, // facesOffsets
        contentType::LEGACY,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateShapesLegacy
(
    const UList<cellShape>& shapes,
    UList<uint8_t>& cellTypes,
    labelUList& vertLabels,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        shapes,
        *this,
        cellTypes,
        vertLabels,
        unused, // offsets
        unused, // faces
        unused, // facesOffsets
        contentType::LEGACY,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateXml
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    labelUList& connectivity,
    labelUList& offsets,
    labelUList& faces,
    labelUList& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::XML,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateShapesXml
(
    const UList<cellShape>& shapes,
    UList<uint8_t>& cellTypes,
    labelUList& connectivity,
    labelUList& offsets,
    labelUList& faces,
    labelUList& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        shapes,
        *this,
        cellTypes,
        connectivity,
        offsets,
        unused, // faces
        unused, // facesOffsets
        contentType::XML,
        maps.cellMap(),
        maps.additionalIds()
    );
}


#undef  definePopulateInternalMethod
#define definePopulateInternalMethod(Type)                                   \
                                                                             \
    void Foam::vtk::vtuSizing::populateInternal                              \
    (                                                                        \
        const polyMesh& mesh,                                                \
        UList<uint8_t>& cellTypes,                                           \
        UList<Type>& connectivity,                                           \
        UList<Type>& offsets,                                                \
        UList<Type>& faces,                                                  \
        UList<Type>& facesOffsets,                                           \
        foamVtkMeshMaps& maps,                                               \
        const enum contentType output                                        \
    ) const                                                                  \
    {                                                                        \
        presizeMaps(maps);                                                   \
                                                                             \
        populateArrays                                                       \
        (                                                                    \
            mesh,                                                            \
            *this,                                                           \
            cellTypes,                                                       \
            connectivity,                                                    \
            offsets,                                                         \
            faces,                                                           \
            facesOffsets,                                                    \
            output,                                                          \
            maps.cellMap(),                                                  \
            maps.additionalIds()                                             \
        );                                                                   \
    }                                                                        \
                                                                             \
    void Foam::vtk::vtuSizing::populateInternal                              \
    (                                                                        \
        const polyMesh& mesh,                                                \
        UList<uint8_t>& cellTypes,                                           \
        UList<Type>& connectivity,                                           \
        UList<Type>& offsets,                                                \
        UList<Type>& faces,                                                  \
        UList<Type>& facesOffsets,                                           \
        labelUList& cellMap,                                                 \
        labelUList& addPointsIds,                                            \
        const enum contentType output                                        \
    ) const                                                                  \
    {                                                                        \
        populateArrays                                                       \
        (                                                                    \
            mesh,                                                            \
            *this,                                                           \
            cellTypes,                                                       \
            connectivity,                                                    \
            offsets,                                                         \
            faces,                                                           \
            facesOffsets,                                                    \
            output,                                                          \
            cellMap,                                                         \
            addPointsIds                                                     \
        );                                                                   \
    }


definePopulateInternalMethod(int);
definePopulateInternalMethod(long);
definePopulateInternalMethod(long long);


#undef definePopulateInternalMethod


// * * * * * * * * * * * * * * Renumber vertices * * * * * * * * * * * * * * //

Foam::labelList Foam::vtk::vtuSizing::copyVertLabelsLegacy
(
    const labelUList& vertLabels,
    const label globalPointOffset
)
{
    if (!globalPointOffset)
    {
        return vertLabels;
    }

    labelList output(vertLabels);
    renumberVertLabelsLegacy(output, globalPointOffset);

    return output;
}


void Foam::vtk::vtuSizing::renumberVertLabelsLegacy
(
    labelUList& vertLabels,
    const label globalPointOffset
)
{
    if (!globalPointOffset)
    {
        return;
    }

    // LEGACY vertLabels = "cells" contains
    // - connectivity
    // [nLabels, vertex labels...]
    // - face-stream
    // [nLabels nFaces, nFace0Pts, id1,id2,..., nFace1Pts, id1,id2,...]

    // Note the simplest volume cell is a tet (4 points, 4 faces)
    // As a poly-face stream this would have
    // 2 for nLabels, nFaces
    // 4 labels (size + ids) per face * 4 == 16 labels
    //
    // Therefore anything with 18 labels or more must be a poly

    auto iter = vertLabels.begin();
    const auto last = vertLabels.end();

    while (iter < last)
    {
        label nLabels = *iter;  // nLabels (for this cell)
        ++iter;

        if (nLabels < 18)
        {
            // Normal primitive type

            while (nLabels--)
            {
                *iter += globalPointOffset;
                ++iter;
            }
        }
        else
        {
            // Polyhedral face-stream (explained above)

            label nFaces = *iter;
            ++iter;

            while (nFaces--)
            {
                nLabels = *iter;  // nLabels (for this face)
                ++iter;

                while (nLabels--)
                {
                    *iter += globalPointOffset;
                    ++iter;
                }
            }
        }
    }
}


Foam::labelList Foam::vtk::vtuSizing::copyVertLabelsXml
(
    const labelUList& vertLabels,
    const label globalPointOffset
)
{
    if (!globalPointOffset)
    {
        return vertLabels;
    }

    labelList output(vertLabels);
    renumberVertLabelsXml(output, globalPointOffset);

    return output;
}


void Foam::vtk::vtuSizing::renumberVertLabelsXml
(
    labelUList& vertLabels,
    const label globalPointOffset
)
{
    if (!globalPointOffset)
    {
        return;
    }

    // XML vertLabels = "connectivity" contains
    // [cell1-verts, cell2-verts, ...]

    for (label& vertId : vertLabels)
    {
        vertId += globalPointOffset;
    }
}


Foam::labelList Foam::vtk::vtuSizing::copyFaceLabelsXml
(
    const labelUList& faceLabels,
    const label globalPointOffset
)
{
    if (!globalPointOffset)
    {
        return faceLabels;
    }

    labelList output(faceLabels);
    renumberFaceLabelsXml(output, globalPointOffset);

    return output;
}


void Foam::vtk::vtuSizing::renumberFaceLabelsXml
(
    labelUList& faceLabels,
    const label globalPointOffset
)
{
    if (!globalPointOffset)
    {
        return;
    }

    // XML face-stream
    // [nFaces, nFace0Pts, id1,id2,..., nFace1Pts, id1,id2,...]

    auto iter = faceLabels.begin();
    const auto last = faceLabels.end();

    while (iter < last)
    {
        label nFaces = *iter;
        ++iter;

        while (nFaces--)
        {
            label nLabels = *iter;
            ++iter;

            while (nLabels--)
            {
                *iter += globalPointOffset;
                ++iter;
            }
        }
    }
}


Foam::labelList Foam::vtk::vtuSizing::copyFaceOffsetsXml
(
    const labelUList& faceOffsets,
    const label prevOffset
)
{
    if (!prevOffset)
    {
        return faceOffsets;
    }

    labelList output(faceOffsets);
    renumberFaceOffsetsXml(output, prevOffset);

    return output;
}


void Foam::vtk::vtuSizing::renumberFaceOffsetsXml
(
    labelUList& faceOffsets,
    const label prevOffset
)
{
    if (!prevOffset)
    {
        return;
    }

    // offsets
    // [-1, off1, off2, ... -1, ..]

    for (label& val : faceOffsets)
    {
        if (val != -1)
        {
            val += prevOffset;
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

void Foam::vtk::vtuSizing::info(Ostream& os) const
{
    os  << "nFieldCells:" << nFieldCells();
    if (nAddCells_)
    {
        os  << " (" << nCells_ << "+" << nAddCells_ << ")";
    }
    else
    {
        os  << " (poly:" << nCellsPoly_ << ")";
    }

    os  << " nFieldPoints:" << nFieldPoints();
    if (nAddPoints_)
    {
        os  << " (" << nPoints_ << "+" << nAddPoints_ << ")";
    }

    os  << " nVertLabels:" << (nVertLabels_ + nAddVerts_);
    if (nAddVerts_)
    {
        os  << " (" << nVertLabels_ << "+" << nAddVerts_ << ")";
    }
    else if (nVertPoly_)
    {
        os  << " (poly:" << nVertPoly_ << ")";
    }

    os << " nFaceLabels:" << nFaceLabels_;
    os << " legacy-count:" << sizeLegacy();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::vtk::vtuSizing::operator==(const vtuSizing& rhs) const
{
    return
    (
        decompose()   == rhs.decompose()
        // required?  && pointOffset() == rhs.pointOffset()
     && nCells()      == rhs.nCells()
     && nPoints()     == rhs.nPoints()
     && nVertLabels() == rhs.nVertLabels()
     && nFaceLabels() == rhs.nFaceLabels()
     && nCellsPoly()  == rhs.nCellsPoly()
     && nVertPoly()   == rhs.nVertPoly()
     && nAddCells()   == rhs.nAddCells()
     && nAddPoints()  == rhs.nAddPoints()
     && nAddVerts()   == rhs.nAddVerts()
    );
}


bool Foam::vtk::vtuSizing::operator!=(const vtuSizing& rhs) const
{
    return !operator==(rhs);
}


// ************************************************************************* //
