/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "rawCoordSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "stringOps.H"
#include "coordSetWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{
    defineTypeName(rawWriter);
    addToRunTimeSelectionTable(coordSetWriter, rawWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, rawWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "rawCoordSetWriterImpl.C"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::rawWriter::rawWriter()
:
    coordSetWriter(),
    streamOpt_(),
    precision_(IOstream::defaultPrecision())
{
    buffering_ = true;
}


Foam::coordSetWriters::rawWriter::rawWriter(const dictionary& options)
:
    coordSetWriter(options),
    streamOpt_
    (
        IOstreamOption::ASCII,
        IOstreamOption::compressionEnum("compression", options)
    ),
    precision_
    (
        options.getOrDefault("precision", IOstream::defaultPrecision())
    )
{
    buffering_ = options.getOrDefault("buffer", true);
}


Foam::coordSetWriters::rawWriter::rawWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    rawWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::rawWriter::rawWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    rawWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::rawWriter::~rawWriter()
{
    close();
}


// * * * * * * * * * * * * * * * * * Controls  * * * * * * * * * * * * * * * //

bool Foam::coordSetWriters::rawWriter::buffering(const bool on)
{
    const bool old(buffering_);
    buffering_ = on;
    return old;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::coordSetWriters::rawWriter::path() const
{
    // Assume !useTracks_, otherwise too fragile

    // 1) rootdir/<TIME>/setName.raw
    // 2) rootdir/setName.raw

    return getExpectedPath("xy");   // Traditionally 'xy', not 'raw'
}


bool Foam::coordSetWriters::rawWriter::writeBuffered()
{
    if (coords_.empty())
    {
        clearBuffers();
        return false;
    }
    const auto& coords = coords_[0];

    // Field:
    // 1) rootdir/<TIME>/setName.raw
    // 2) rootdir/setName.raw

    fileName outputFile = path();

    if (!isDir(outputFile.path()))
    {
        mkDir(outputFile.path());
    }

    OFstream os(outputFile, streamOpt_);
    os.precision(precision_);

    writeBufferContents(os, coords, " \t");

    clearBuffers();

    return true;
}


template<class Type>
Foam::fileName Foam::coordSetWriters::rawWriter::writeTemplate
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

    if (useTracks_ || !buffering_)
    {
        UPtrList<const Field<Type>> fieldPtrs(repackageFields(values));
        return writeTemplate(fieldName, fieldPtrs);
    }

    // Buffering version
    appendField(fieldName, values);
    return path();
}


template<class Type>
Foam::fileName Foam::coordSetWriters::rawWriter::writeTemplate
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
defineCoordSetWriterWriteFields(Foam::coordSetWriters::rawWriter);


// ************************************************************************* //
