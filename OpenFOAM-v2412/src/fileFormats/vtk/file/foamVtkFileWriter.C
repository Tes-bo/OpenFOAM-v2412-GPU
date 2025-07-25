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

#include "foamVtkFileWriter.H"
#include "globalIndex.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::vtk::fileWriter::outputState
>
Foam::vtk::fileWriter::stateNames
({
    { outputState::CLOSED, "closed" },
    { outputState::OPENED, "opened" },
    { outputState::DECLARED, "declared" },
    { outputState::FIELD_DATA, "FieldData" },
    { outputState::PIECE, "Piece" },
    { outputState::CELL_DATA, "CellData" },
    { outputState::POINT_DATA, "PointData" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::vtk::fileWriter::checkFormatterValidity() const
{
    // In parallel can be unallocated on non-master nodes
    if ((parallel_ ? Pstream::master() : true) && !format_)
    {
        FatalErrorInFunction
            << "unallocated formatter" << endl
            << exit(FatalError);
    }
}


Foam::Ostream& Foam::vtk::fileWriter::reportBadState
(
    Ostream& os,
    outputState expected
) const
{
    os  << "Bad writer state (" << stateNames[state_]
        << ") - should be (" << stateNames[expected] << ')';
    return os;
}


Foam::Ostream& Foam::vtk::fileWriter::reportBadState
(
    Ostream& os,
    outputState expected1,
    outputState expected2
) const
{
    reportBadState(os, expected1)
        << " or (" << stateNames[expected2] << ')';
    return os;
}


bool Foam::vtk::fileWriter::enter_Piece()
{
    // Finish other output
    endFieldData();

    if (isState(outputState::OPENED))
    {
        beginFile();
    }
    if (notState(outputState::DECLARED))
    {
        reportBadState(FatalErrorInFunction, outputState::DECLARED)
            << exit(FatalError);
    }
    state_ = outputState::PIECE;
    nCellData_ = nPointData_ = 0;

    return true;
}


bool Foam::vtk::fileWriter::endPiece()
{
    // Finish other output
    endCellData();
    endPointData();

    if (notState(outputState::PIECE))
    {
        // Skip if not in Piece
        return false;
    }
    state_ = outputState::DECLARED; // Mark as having been flushed

    if (format_)
    {
        format().endPiece();
    }

    return true;
}


bool Foam::vtk::fileWriter::enter_CellData(label nEntries, label nFields)
{
    // Already in CellData?
    if (isState(outputState::CELL_DATA)) return false;

    // Finish other output
    endPointData();

    if (notState(outputState::PIECE))
    {
        reportBadState(FatalErrorInFunction, outputState::PIECE)
            << exit(FatalError);
    }

    nCellData_ = 0;

    // Do nothing for legacy when nFields == 0
    if (legacy() && !nFields) return false;

    state_ = outputState::CELL_DATA;

    if (format_)
    {
        if (legacy())
        {
            legacy::beginCellData(format(), nEntries, nFields);
        }
        else
        {
            format().beginCellData();
        }
    }

    return true;
}


bool Foam::vtk::fileWriter::enter_PointData(label nEntries, label nFields)
{
    // Already in PointData?
    if (isState(outputState::POINT_DATA)) return false;

    // Finish other output
    endCellData();

    if (notState(outputState::PIECE))
    {
        reportBadState(FatalErrorInFunction, outputState::PIECE)
            << exit(FatalError);
    }

    nPointData_ = 0;

    // Do nothing for legacy when nFields == 0
    if (legacy() && !nFields) return false;

    state_ = outputState::POINT_DATA;

    if (format_)
    {
        if (legacy())
        {
            legacy::beginPointData(format(), nEntries, nFields);
        }
        else
        {
            format().beginPointData();
        }
    }

    return true;
}


void Foam::vtk::fileWriter::endDataArray()
{
    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


void Foam::vtk::fileWriter::beginPoints(const label nPoints)
{
    if (format_)
    {
        if (legacy())
        {
            legacy::beginPoints(os_, nPoints);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, 3>(nPoints);

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float, 3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);
        }
    }
}


void Foam::vtk::fileWriter::endPoints()
{
    if (format_)
    {
        format().flush();
        format().endDataArray();

        if (!legacy())
        {
            format()
                .endTag(vtk::fileTag::POINTS);
        }
    }
}


bool Foam::vtk::fileWriter::exit_File()
{
    // Finish other output
    endFieldData();
    endPiece();

    if (isState(outputState::DECLARED))
    {
        if (format_ && !legacy())
        {
            format().endTag(contentType_).endVTKFile();
        }
        state_ = outputState::OPENED;  // Mark as having been flushed
    }

    // Must now be in CLOSED or OPENED states only

    if (isState(outputState::CLOSED) || isState(outputState::OPENED))
    {
        return true;
    }

    reportBadState(WarningInFunction, outputState::CLOSED, outputState::OPENED)
        << " for contentType (" << vtk::fileTagNames[contentType_] << ')'
        << nl << endl;

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::fileWriter::fileWriter
(
    const vtk::fileTag contentType,
    const vtk::outputOptions opts
)
:
    state_(outputState::CLOSED),
    contentType_(contentType),
    parallel_(false),
    opts_(opts),
    nCellData_(0),
    nPointData_(0),
    outputFile_(),
    format_(nullptr),
    os_()
{
    // We do not currently support append mode at all
    opts_.append(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::fileWriter::~fileWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::fileWriter::open(const fileName& file, bool parallel)
{
    if (notState(outputState::CLOSED))
    {
        reportBadState(FatalErrorInFunction, outputState::CLOSED)
            << exit(FatalError);
    }

    if (format_)
    {
        format_.reset(nullptr);
        os_.close();
    }
    nCellData_ = nPointData_ = 0;
    outputFile_ = file;

    if
    (
        legacy()
      ? outputFile_.has_ext(vtk::fileExtension[contentType_])
      : outputFile_.has_ext(vtk::legacy::fileExtension)
    )
    {
        // Inappropriate extension. Legacy instead of xml, or vice versa.
        outputFile_.remove_ext();
    }

    if (!outputFile_.has_ext(ext()))
    {
        // Add extension if required
        outputFile_.ext(ext());
    }


    // Only set parallel flag if really is a parallel run.
    parallel_ = parallel && Pstream::parRun();

    // Open a file and attach a formatter
    // - on master (always)
    // - on subproc (if not parallel)
    //
    // This means we can always check if format_ is defined to know if output
    // is desired on any particular process.

    if (!parallel_ || Pstream::master())
    {
        mkDir(outputFile_.path());

        os_.open(outputFile_);

        format_ = opts_.newFormatter(os_);
    }

    state_ = outputState::OPENED;
    return true;
}


void Foam::vtk::fileWriter::close()
{
    exit_File();

    if (format_)
    {
        format_.reset(nullptr);
        os_.close();
    }

    state_ = outputState::CLOSED;
    outputFile_.clear();
    nCellData_ = nPointData_ = 0;
}


bool Foam::vtk::fileWriter::beginFile(std::string title)
{
    if (isState(outputState::DECLARED))
    {
        // Skip if already emitted
        return false;
    }
    if (notState(outputState::OPENED))
    {
        reportBadState(FatalErrorInFunction, outputState::OPENED)
            << exit(FatalError);
    }
    state_ = outputState::DECLARED;

    if (format_)
    {
        if (legacy())
        {
            legacy::fileHeader(format(), title, contentType_);
        }
        else
        {
            // XML (inline)

            format().xmlHeader();

            if (title.size())
            {
                format().xmlComment(title);
            }

            format().beginVTKFile(contentType_);
        }
    }

    return true;
}


bool Foam::vtk::fileWriter::beginFieldData(label nFields)
{
    // Do nothing for legacy when nFields == 0
    if (legacy() && !nFields) return false;

    if (isState(outputState::OPENED))
    {
        beginFile();
    }
    if (notState(outputState::DECLARED))
    {
        reportBadState(FatalErrorInFunction, outputState::DECLARED)
            << exit(FatalError);
    }
    state_ = outputState::FIELD_DATA;

    if (format_)
    {
        if (legacy())
        {
            legacy::beginFieldData(format(), nFields);
        }
        else
        {
            format().beginFieldData();
        }
    }

    return true;
}


bool Foam::vtk::fileWriter::endFieldData()
{
    if (notState(outputState::FIELD_DATA))
    {
        // Skip if not in FieldData
        return false;
    }
    state_ = outputState::DECLARED; // Toggle back to DECLARED

    if (format_ && !legacy())
    {
        format().endFieldData();
    }

    return true;
}


bool Foam::vtk::fileWriter::endCellData()
{
    if (notState(outputState::CELL_DATA))
    {
        // Skip if not in CellData
        return false;
    }
    state_ = outputState::PIECE; // Toggle back to PIECE

    if (format_ && !legacy())
    {
        format().endCellData();
    }

    return true;
}


bool Foam::vtk::fileWriter::endPointData()
{
    if (notState(outputState::POINT_DATA))
    {
        // Skip if not in PointData
        return false;
    }
    state_ = outputState::PIECE; // Toggle back to PIECE

    if (format_ && !legacy())
    {
        format().endPointData();
    }

    return true;
}


void Foam::vtk::fileWriter::writeTimeValue(scalar timeValue)
{
    // Convenience - switch to FieldData
    if (isState(outputState::OPENED) || isState(outputState::DECLARED))
    {
        beginFieldData(1);
    }
    if (notState(outputState::FIELD_DATA))
    {
        reportBadState(FatalErrorInFunction, outputState::FIELD_DATA)
            << exit(FatalError);
    }

    // No collectives - can skip on sub-procs
    if (!format_) return;

    if (legacy())
    {
        legacy::writeTimeValue(format(), timeValue);
    }
    else
    {
        format().writeTimeValue(timeValue);
    }
}


bool Foam::vtk::fileWriter::writeProcIDs(const label nValues)
{
    // Write procIDs whenever running in parallel

    if (!Pstream::parRun())
    {
        return false;  // Non-parallel: skip
    }

    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState
        (
            FatalErrorInFunction,
            outputState::CELL_DATA,
            outputState::POINT_DATA
        )   << " for procID field" << nl << endl
            << exit(FatalError);

        return false;
    }


    const globalIndex procAddr
    (
        parallel_
      ? globalIndex(globalIndex::gatherOnly{}, nValues)
      : globalIndex(globalIndex::gatherNone{}, nValues)
    );

    const label totalCount = procAddr.totalSize();

    this->beginDataArray<label>("procID", totalCount);

    bool good = false;

    if (parallel_)
    {
        if (Pstream::master())
        {
            // Per-processor ids
            for (const label proci : procAddr.allProcs())
            {
                vtk::write(format(), proci, procAddr.localSize(proci));
            }
            good = true;
        }
    }
    else
    {
        vtk::write(format(), label(Pstream::myProcNo()), totalCount);
        good = true;
    }


    this->endDataArray();

    // MPI barrier
    return parallel_ ? returnReduceOr(good) : good;
}


// ************************************************************************* //
