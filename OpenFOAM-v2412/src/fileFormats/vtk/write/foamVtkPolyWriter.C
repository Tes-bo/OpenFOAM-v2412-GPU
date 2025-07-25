/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "foamVtkPolyWriter.H"
#include "foamVtkOutput.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// The connectivity count for a list of edges
static inline label countConnectivity(const edgeList& edges)
{
    return 2 * edges.size();  // An edge always has two ends
}


// The connectivity count for a list of faces
static label countConnectivity(const faceList& faces)
{
    label nConnectivity = 0;

    for (const face& f : faces)
    {
        nConnectivity += f.size();
    }

    return nConnectivity;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::polyWriter::beginPiece
(
    const pointField& points,
    const edgeList& edges
)
{
    // Basic sizes
    nLocalPoints_ = points.size();
    nLocalVerts_  = 0;
    nLocalLines_  = edges.size();
    nLocalPolys_  = 0;

    numberOfPoints_ = nLocalPoints_;
    numberOfCells_  = nLocalLines_;

    if (parallel_)
    {
        reduce(numberOfPoints_, sumOp<label>());
        reduce(numberOfCells_,  sumOp<label>());
    }


    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        format().tag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, numberOfPoints_,
            vtk::fileAttr::NUMBER_OF_LINES,  numberOfCells_
            // AND: vtk::fileAttr::NUMBER_OF_POLYS,  0
        );
    }
}


void Foam::vtk::polyWriter::beginPiece
(
    const pointField& points,
    const faceList& faces
)
{
    // Basic sizes
    nLocalPoints_ = points.size();
    nLocalVerts_  = 0;
    nLocalLines_  = 0;
    nLocalPolys_  = faces.size();

    numberOfPoints_ = nLocalPoints_;
    numberOfCells_  = nLocalPolys_;

    if (parallel_)
    {
        reduce(numberOfPoints_, sumOp<label>());
        reduce(numberOfCells_,  sumOp<label>());
    }


    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        format().tag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, numberOfPoints_,
            vtk::fileAttr::NUMBER_OF_POLYS,  numberOfCells_
            // AND: vtk::fileAttr::NUMBER_OF_LINES,  0
        );
    }
}


void Foam::vtk::polyWriter::writePoints
(
    const pointField& points
)
{
    this->beginPoints(numberOfPoints_);

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), points);
    }
    else
    {
        vtk::writeList(format(), points);

    }

    this->endPoints();
}


void Foam::vtk::polyWriter::writeLinesLegacy
(
    const edgeList& edges,
    const label pointOffset
)
{
    // Connectivity count without additional storage (done internally)
    const label nLocalConns = countConnectivity(edges);

    label nLines = nLocalLines_;
    label nConns = nLocalConns;

    if (parallel_)
    {
        reduce(nLines, sumOp<label>());
        reduce(nConns, sumOp<label>());
    }

    if (nLines != numberOfCells_)
    {
        FatalErrorInFunction
            << "Expecting " << numberOfCells_
            << " edges, but found " << nLines
            << exit(FatalError);
    }

    legacy::beginLines(os_, nLines, nConns);

    labelList vertLabels(nLocalLines_ + nLocalConns);

    {
        // Legacy: size + connectivity together
        // [nPts, id1, id2, ..., nPts, id1, id2, ...]

        auto iter = vertLabels.begin();

        const label off = pointOffset;

        for (const edge& e : edges)
        {
            *iter = e.size();   // The size prefix (always 2 for an edge)
            ++iter;

            *iter = off + e.first();    // Vertex labels
            ++iter;

            *iter = off + e.second();
            ++iter;
        }
    }


    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), vertLabels);
    }
    else
    {
        vtk::writeList(format(), vertLabels);
    }

    if (format_)
    {
        format().flush();
    }
}


void Foam::vtk::polyWriter::writeLines
(
    const edgeList& edges,
    const label pointOffset
)
{
    // Connectivity count without additional storage (done internally)
    const label nLocalConns = countConnectivity(edges);

    if (format_)
    {
        format().tag(vtk::fileTag::LINES);
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalConns);

        label nConns = nLocalConns;

        if (parallel_)
        {
            reduce(nConns, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nConns);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad * sizeof(label));
        }

        {
            // XML: connectivity only
            // [id1, id2, ..., id1, id2, ...]

            auto iter = vertLabels.begin();

            const label off = pointOffset;

            for (const edge& e : edges)
            {
                // Edge vertex labels
                *iter = off + e.first();
                ++iter;

                *iter = off + e.second();
                ++iter;
            }
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertLabels);
        }
        else
        {
            vtk::writeList(format(), vertLabels);
        }

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        labelList vertOffsets(nLocalLines_);
        label nOffs = vertOffsets.size();

        if (parallel_)
        {
            reduce(nOffs, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }


        // processor-local connectivity offsets
        label off =
        (
            parallel_ ? globalIndex::calcOffset(nLocalConns) : 0
        );


        auto iter = vertOffsets.begin();

        for (const edge& e : edges)
        {
            off += e.size();   // End offset
            *iter = off;
            ++iter;
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertOffsets);
        }
        else
        {
            vtk::writeList(format_.ref(), vertOffsets);
        }


        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    if (format_)
    {
        format().endTag(vtk::fileTag::LINES);
    }
}


void Foam::vtk::polyWriter::writePolysLegacy
(
    const faceList& faces,
    const label pointOffset
)
{
    // Connectivity count without additional storage (done internally)
    const label nLocalConns = countConnectivity(faces);

    label nPolys = nLocalPolys_;
    label nConns = nLocalConns;

    if (parallel_)
    {
        reduce(nPolys, sumOp<label>());
        reduce(nConns, sumOp<label>());
    }

    if (nPolys != numberOfCells_)
    {
        FatalErrorInFunction
            << "Expecting " << numberOfCells_
            << " faces, but found " << nPolys
            << exit(FatalError);
    }

    legacy::beginPolys(os_, nPolys, nConns);

    labelList vertLabels(nLocalPolys_ + nLocalConns);

    {
        // Legacy: size + connectivity together
        // [nPts, id1, id2, ..., nPts, id1, id2, ...]

        auto iter = vertLabels.begin();

        const label off = pointOffset;

        for (const face& f : faces)
        {
            *iter = f.size();       // The size prefix
            ++iter;

            for (const label id : f)
            {
                *iter = id + off;   // Vertex label
                ++iter;
            }
        }
    }


    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), vertLabels);
    }
    else
    {
        vtk::writeList(format(), vertLabels);
    }

    if (format_)
    {
        format().flush();
    }
}


void Foam::vtk::polyWriter::writePolys
(
    const faceList& faces,
    const label pointOffset
)
{
    // Connectivity count without additional storage (done internally)
    const label nLocalConns = countConnectivity(faces);

    if (format_)
    {
        format().tag(vtk::fileTag::POLYS);
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalConns);

        label nConns = nLocalConns;

        if (parallel_)
        {
            reduce(nConns, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nConns);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad * sizeof(label));
        }

        {
            // XML: connectivity only
            // [id1, id2, ..., id1, id2, ...]

            auto iter = vertLabels.begin();

            label off = pointOffset;

            for (const face& f : faces)
            {
                for (const label id : f)
                {
                    *iter = id + off;  // Face vertex label
                    ++iter;
                }
            }
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertLabels);
        }
        else
        {
            vtk::writeList(format(), vertLabels);
        }

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        labelList vertOffsets(nLocalPolys_);
        label nOffs = vertOffsets.size();

        if (parallel_)
        {
            reduce(nOffs, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }


        // processor-local connectivity offsets
        label off =
        (
            parallel_ ? globalIndex::calcOffset(nLocalConns) : 0
        );


        auto iter = vertOffsets.begin();

        for (const face& f : faces)
        {
            off += f.size();  // End offset
            *iter = off;
            ++iter;
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertOffsets);
        }
        else
        {
            vtk::writeList(format_.ref(), vertOffsets);
        }


        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    if (format_)
    {
        format().endTag(vtk::fileTag::POLYS);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::polyWriter::polyWriter
(
    const vtk::outputOptions opts
)
:
    vtk::fileWriter(vtk::fileTag::POLY_DATA, opts),
    numberOfPoints_(0),
    numberOfCells_(0),
    nLocalPoints_(0),
    nLocalVerts_(0),
    nLocalLines_(0),
    nLocalPolys_(0)
{
    // We do not currently support append mode
    opts_.append(false);
}


Foam::vtk::polyWriter::polyWriter
(
    const fileName& file,
    bool parallel
)
:
    // Default parameter fails for gcc-4.8.5, thus specify format here
    polyWriter(vtk::formatType::INLINE_BASE64)
{
    open(file, parallel);
}


Foam::vtk::polyWriter::polyWriter
(
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    polyWriter(opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::polyWriter::writeGeometry()
{
    FatalErrorInFunction
        << "Method was not overloaded, called without a geometry!!" << nl
        << "    Indicates a programming error" << nl << endl
        << abort(FatalError);

    return false;
}


bool Foam::vtk::polyWriter::writeLineGeometry
(
    const pointField& points,
    const edgeList& edges
)
{
    enter_Piece();

    beginPiece(points, edges);

    writePoints(points);

    const label pointOffset =
    (
        parallel_ ? globalIndex::calcOffset(nLocalPoints_) : 0
    );

    if (legacy())
    {
        writeLinesLegacy(edges, pointOffset);
    }
    else
    {
        writeLines(edges, pointOffset);
    }

    return true;
}


bool Foam::vtk::polyWriter::writePolyGeometry
(
    const pointField& points,
    const faceList& faces
)
{
    enter_Piece();

    beginPiece(points, faces);

    writePoints(points);

    const label pointOffset =
    (
        parallel_ ? globalIndex::calcOffset(nLocalPoints_) : 0
    );

    if (legacy())
    {
        writePolysLegacy(faces, pointOffset);
    }
    else
    {
        writePolys(faces, pointOffset);
    }

    return true;
}


bool Foam::vtk::polyWriter::beginCellData(label nFields)
{
    return enter_CellData(numberOfCells_, nFields);
}


bool Foam::vtk::polyWriter::beginPointData(label nFields)
{
    return enter_PointData(numberOfPoints_, nFields);
}


// ************************************************************************* //
