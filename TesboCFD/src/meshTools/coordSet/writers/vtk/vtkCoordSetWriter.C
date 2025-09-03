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

#include "vtkCoordSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "foamVtkCoordSetWriter.H"
#include "coordSetWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{
    defineTypeName(vtkWriter);
    addToRunTimeSelectionTable(coordSetWriter, vtkWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, vtkWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::vtkWriter::vtkWriter()
:
    coordSetWriter(),
    fmtType_(static_cast<unsigned>(vtk::formatType::INLINE_BASE64)),
    precision_(IOstream::defaultPrecision()),
    writer_(nullptr)
{}


Foam::coordSetWriters::vtkWriter::vtkWriter
(
    const vtk::outputOptions& opts
)
:
    coordSetWriter(),
    fmtType_(static_cast<unsigned>(opts.fmt())),
    precision_(opts.precision()),
    writer_(nullptr)
{}


Foam::coordSetWriters::vtkWriter::vtkWriter(const dictionary& options)
:
    coordSetWriter(options),
    fmtType_(static_cast<unsigned>(vtk::formatType::INLINE_BASE64)),
    precision_
    (
        options.getOrDefault("precision", IOstream::defaultPrecision())
    ),
    writer_(nullptr)
{
    // format: ascii | binary
    // legacy: true | false

    vtk::outputOptions opts(vtk::formatType::INLINE_BASE64);
    opts.ascii
    (
        IOstreamOption::ASCII
     == IOstreamOption::formatEnum("format", options, IOstreamOption::BINARY)
    );

    opts.legacy(options.getOrDefault("legacy", false));

    // Convert back to raw data type
    fmtType_ = static_cast<unsigned>(opts.fmt());
}


Foam::coordSetWriters::vtkWriter::vtkWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    vtkWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::vtkWriter::vtkWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    vtkWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::vtkWriter::~vtkWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::coordSetWriters::vtkWriter::path() const
{
    // 1) rootdir/<TIME>/setName.{vtk|vtp}
    // 2) rootdir/setName.{vtk|vtp}

    // From raw unsigned value to vtk::outputOptions
    vtk::outputOptions opts(static_cast<vtk::formatType>(fmtType_));

    return getExpectedPath(vtk::polyWriter::ext(opts));
}


void Foam::coordSetWriters::vtkWriter::close(bool force)
{
    writer_.clear();
    coordSetWriter::close(force);
}


void Foam::coordSetWriters::vtkWriter::beginTime(const Time& t)
{
    writer_.clear();
    coordSetWriter::beginTime(t);
}


void Foam::coordSetWriters::vtkWriter::beginTime(const instant& inst)
{
    writer_.clear();
    coordSetWriter::beginTime(inst);
}


void Foam::coordSetWriters::vtkWriter::endTime()
{
    writer_.clear();
    coordSetWriter::endTime();
}


Foam::fileName Foam::coordSetWriters::vtkWriter::write()
{
    checkOpen();
    if (needsUpdate())
    {
        writer_.clear();
    }
    merge();

    if (coords_.empty())
    {
        return fileName::null;
    }

    // From raw unsigned values to vtk::outputOptions
    vtk::outputOptions opts(static_cast<vtk::formatType>(fmtType_), precision_);


    // Geometry:  rootdir/<TIME>/setName.{vtk|vtp}

    fileName outputFile = getExpectedPath(vtk::polyWriter::ext(opts));

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    if (!writer_ && true)  // always (non-parallel)
    {
        UPtrList<const pointField> points(coords_.size());
        forAll(coords_, tracki)
        {
            points.set(tracki, coords_.get(tracki));
        }

        writer_.reset
        (
            new vtk::coordSetWriter
            (
                points,
                opts,
                outputFile,
                false  // serial!
            )
        );

        if (useTracks_ || coords_.size() > 1)
        {
            writer_->setElementType(vtk::coordSetWriter::LINE_ELEMENTS);
        }

        if (this->hasTime())
        {
            // Time name in title
            writer_->setTime(currTime_);
            writer_->writeTimeValue();
        }
        else
        {
            // Set name in title
            writer_->beginFile(outputPath_.stem());
        }

        writer_->writeGeometry();
    }

    wroteGeom_ = true;
    return outputFile;
}


// * * * * * * * * * * * * * * * Implementation * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::coordSetWriters::vtkWriter::writeTemplate
(
    const word& fieldName,
    const UPtrList<const Field<Type>>& fieldPtrs
)
{
    if (coords_.size() != fieldPtrs.size())
    {
        FatalErrorInFunction
            << "Attempted to write field: " << fieldName
            << " (" << fieldPtrs.size() << " entries) for "
            << coords_.size() << " sets" << nl
            << exit(FatalError);
    }

    // Open file, writing geometry (if required)
    fileName outputFile = this->write();

    if (!nFields_ && writer_->legacy())
    {
        // Emit error message, but attempt to recover anyhow
        nFields_ = 1;

        FatalErrorInFunction
            << "Using VTK legacy format, but did not define nFields!"
            << nl
            << "Assuming nFields=1 (may be incorrect) and continuing..."
            << nl
            << "    Field " << fieldName << " to " << outputFile << nl;

        Info<< FatalError;
        Info<< endl;
    }

    writer_->beginPointData(nFields_);

    writer_->writePointData(fieldName, fieldPtrs);

    wroteGeom_ = true;
    return outputFile;
}


template<class Type>
Foam::fileName Foam::coordSetWriters::vtkWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& values
)
{
    checkOpen();
    if (coords_.empty())
    {
        return fileName::null;
    }

    UPtrList<const Field<Type>> fieldPtrs(repackageFields(values));
    return writeTemplate(fieldName, fieldPtrs);
}


template<class Type>
Foam::fileName Foam::coordSetWriters::vtkWriter::writeTemplate
(
    const word& fieldName,
    const List<Field<Type>>& fieldValues
)
{
    checkOpen();
    if (coords_.empty())
    {
        return fileName::null;
    }
    useTracks_ = true;  // Extra safety

    UPtrList<const Field<Type>> fieldPtrs(repackageFields(fieldValues));
    return writeTemplate(fieldName, fieldPtrs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineCoordSetWriterWriteFields(Foam::coordSetWriters::vtkWriter);


// ************************************************************************* //
