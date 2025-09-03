/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "foamVtkCoordSetWriter.H"
#include "foamVtkOutput.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::coordSetWriter::beginPiece()
{
    // Basic sizes
    nLocalPoints_ = 0;
    nLocalVerts_  = 0;
    nLocalLines_  = 0;
    nLocalPolys_  = 0;

    for (const pointField& pts : points_)
    {
        const label npts = pts.size();
        nLocalPoints_ += npts;

        if (npts)
        {
            ++nLocalLines_;
        }
    }

    switch (elemOutput_)
    {
        case elemOutputType::NO_ELEMENTS:
        {
            nLocalVerts_ = nLocalLines_ = 0;
            break;
        }
        case elemOutputType::DEFAULT_ELEMENTS:
        {
            if (points_.size() < 2)
            {
                //OR  nLocalVerts_ = nLocalPoints_;
                nLocalVerts_ = 0;
                nLocalLines_ = 0;
            }
            break;
        }
        case elemOutputType::POINT_ELEMENTS:
        {
            nLocalVerts_ = nLocalPoints_;
            nLocalLines_ = 0;
            break;
        }
        case elemOutputType::LINE_ELEMENTS:
        {
            // Already determined
            break;
        }
    }

    // Update sizes, similar to
    // vtk::polyWriter::beginPiece(const pointField&, const edgeList&)

    numberOfPoints_ = nLocalPoints_;
    numberOfCells_  = nLocalLines_;

    // if (parallel_)
    // {
    //     reduce(numberOfPoints_, sumOp<label>());
    //     reduce(numberOfCells_,  sumOp<label>());
    // }


    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        format().openTag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, numberOfPoints_
        );
        if (nLocalVerts_)
        {
            format().xmlAttr(vtk::fileAttr::NUMBER_OF_VERTS, nLocalVerts_);
        }
        if (nLocalLines_)
        {
            format().xmlAttr(vtk::fileAttr::NUMBER_OF_LINES, nLocalLines_);
        }
        format().closeTag();
    }
}


void Foam::vtk::coordSetWriter::writePoints()
{
    this->beginPoints(numberOfPoints_);  //<- same as nLocalPoints_

    {
        for (const pointField& pts : points_)
        {
            vtk::writeList(format(), pts);
        }
    }

    this->endPoints();
}


void Foam::vtk::coordSetWriter::writeVertsLegacy()
{
    if (!nLocalVerts_)
    {
        return;  // Nothing to do
    }

    // connectivity = 1 per vertex
    const label nLocalConns = nLocalVerts_;

    legacy::beginVerts(os_, nLocalVerts_, nLocalConns);

    labelList vertLabels(nLocalVerts_ + nLocalConns);

    auto iter = vertLabels.begin();

    for (label pointi = 0; pointi < nLocalVerts_; ++pointi)
    {
        *iter++ = 1;
        *iter++ = pointi;
    }

    vtk::writeList(format(), vertLabels);

    if (format_)
    {
        format().flush();
    }
}


void Foam::vtk::coordSetWriter::writeLinesLegacy()
{
    if (!nLocalLines_)
    {
        return;  // Nothing to do
    }

    // connectivity = use each point
    label nLocalConns = nLocalPoints_;

    legacy::beginLines(os_, nLocalLines_, nLocalConns);

    labelList vertLabels(nLocalLines_ + nLocalConns);

    auto iter = vertLabels.begin();

    label localPointi = 0;
    for (const pointField& pts : points_)
    {
        label npts = pts.size();

        if (npts)
        {
            *iter++ = npts;
            while (npts--)
            {
                *iter++ = localPointi;
                ++localPointi;
            }
        }
    }

    vtk::writeList(format(), vertLabels);

    if (format_)
    {
        format().flush();
    }
}


void Foam::vtk::coordSetWriter::writeVerts()
{
    if (!nLocalVerts_)
    {
        return;  // Nothing to do
    }

    // connectivity = 1 per vertex
    const label nLocalConns = nLocalVerts_;

    if (format_)
    {
        format().tag(vtk::fileTag::VERTS);
    }

    //
    // 'offsets'  (connectivity offsets)
    //
    {
        labelList vertOffsets(nLocalVerts_);
        label nOffs = vertOffsets.size();

        // if (parallel_)
        // {
        //     reduce(nOffs, sumOp<label>());
        // }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }

        // processor-local connectivity offsets

        label off = 0;

        /// label off =
        /// (
        ///     parallel_ ? globalIndex::calcOffset(nLocalConns) : 0
        /// );

        auto iter = vertOffsets.begin();

        for (label pointi = 0; pointi < nLocalVerts_; ++pointi)
        {
            off += 1;  // End offset
            *iter = off;
            ++iter;
        }

        vtk::writeList(format_.ref(), vertOffsets);

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalConns);

        label nConns = nLocalConns;

        // if (parallel_)
        // {
        //     reduce(nConns, sumOp<label>());
        // }

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

            for (label pointi = 0; pointi < nLocalVerts_; ++pointi)
            {
                *iter++ = pointi;
            }
        }

        vtk::writeList(format(), vertLabels);

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }


    if (format_)
    {
        format().endTag(vtk::fileTag::VERTS);
    }
}


void Foam::vtk::coordSetWriter::writeLines()
{
    if (!nLocalLines_)
    {
        return;  // Nothing to do
    }

    // connectivity = use each point
    label nLocalConns = nLocalPoints_;

    if (format_)
    {
        format().tag(vtk::fileTag::LINES);
    }

    //
    // 'offsets'  (connectivity offsets)
    //
    {
        labelList vertOffsets(nLocalLines_);
        label nOffs = vertOffsets.size();

        // if (parallel_)
        // {
        //     reduce(nOffs, sumOp<label>());
        // }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }

        // processor-local connectivity offsets

        label off = 0;

        /// label off =
        /// (
        ///     parallel_ ? globalIndex::calcOffset(nLocalConns) : 0
        /// );

        auto iter = vertOffsets.begin();

        for (const pointField& pts : points_)
        {
            const label npts = pts.size();

            if (npts)
            {
                off += npts;  // End offset
                *iter = off;
                ++iter;
            }
        }

        vtk::writeList(format_.ref(), vertOffsets);

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalConns);

        label nConns = nLocalConns;

        // if (parallel_)
        // {
        //     reduce(nConns, sumOp<label>());
        // }

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

            label localPointi = 0;
            for (const pointField& pts : points_)
            {
                label npts = pts.size();

                while (npts--)
                {
                    *iter++ = localPointi;
                    ++localPointi;
                }
            }
        }


        vtk::writeList(format(), vertLabels);

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::coordSetWriter::coordSetWriter
(
    const UPtrList<const pointField>& points,
    const vtk::outputOptions opts
)
:
    vtk::polyWriter(opts),

    points_(points),
    instant_(),
    elemOutput_(DEFAULT_ELEMENTS)
{}


Foam::vtk::coordSetWriter::coordSetWriter
(
    const UPtrList<const pointField>& points,
    const fileName& file,
    bool parallel
)
:
    coordSetWriter(points)
{
    open(file, parallel);
}


Foam::vtk::coordSetWriter::coordSetWriter
(
    const UPtrList<const pointField>& points,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    coordSetWriter(points, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::coordSetWriter::setElementType
(
    const elemOutputType elemOutput
)
{
    elemOutput_ = elemOutput;
}


bool Foam::vtk::coordSetWriter::open
(
    const fileName& file,
    bool parallel
)
{
    return vtk::polyWriter::open(file, false);  // non-parallel only
}


void Foam::vtk::coordSetWriter::setTime(const instant& inst)
{
    instant_ = inst;
}


bool Foam::vtk::coordSetWriter::beginFile(std::string title)
{
    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }

    if (!instant_.name().empty())
    {
        return vtk::fileWriter::beginFile
        (
            "time='" + instant_.name() + "'"
        );
    }

    // Provide default title
    return vtk::fileWriter::beginFile("coord-set");
}


bool Foam::vtk::coordSetWriter::writeGeometry()
{
    enter_Piece();

    beginPiece();

    writePoints();

    //const label pointOffset =
    //(
    //   parallel_ ? globalIndex::calcOffset(nLocalPoints_) : 0
    //);

    if (legacy())
    {
        writeVertsLegacy();
        writeLinesLegacy();
    }
    else
    {
        writeVerts();
        writeLines();
    }

    return true;
}


void Foam::vtk::coordSetWriter::writeTimeValue()
{
    if (!instant_.name().empty())
    {
        vtk::fileWriter::writeTimeValue(instant_.value());
    }
}


void Foam::vtk::coordSetWriter::piece
(
    const UPtrList<const pointField>& points
)
{
    endPiece();

    points_ = points;
}


bool Foam::vtk::coordSetWriter::writeProcIDs()
{
    // Ignore
    return false;
}


// ************************************************************************* //
