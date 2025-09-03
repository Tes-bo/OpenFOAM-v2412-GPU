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

#include "xmgraceCoordSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "coordSetWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{
    defineTypeName(xmgraceWriter);
    addToRunTimeSelectionTable(coordSetWriter, xmgraceWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, xmgraceWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "xmgraceCoordSetWriterImpl.C"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::xmgraceWriter::xmgraceWriter()
:
    coordSetWriter(),
    streamOpt_(),
    precision_(IOstream::defaultPrecision()),
    ofile_(nullptr),
    nWritten_(0)
{
    buffering_ = true;
}


Foam::coordSetWriters::xmgraceWriter::xmgraceWriter(const dictionary& options)
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
    ),
    ofile_(nullptr),
    nWritten_(0)
{
    buffering_ = options.getOrDefault("buffer", true);
}


Foam::coordSetWriters::xmgraceWriter::xmgraceWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    xmgraceWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::xmgraceWriter::xmgraceWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    xmgraceWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::xmgraceWriter::~xmgraceWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSetWriters::xmgraceWriter::buffering(const bool on)
{
    const bool old(buffering_);
    buffering_ = on;
    return old;
}


Foam::fileName Foam::coordSetWriters::xmgraceWriter::path() const
{
    // Assume !useTracks_, otherwise too fragile

    // 1) rootdir/<TIME>/setName.agr
    // 2) rootdir/setName.agr

    return getExpectedPath("agr");
}


void Foam::coordSetWriters::xmgraceWriter::close(bool force)
{
    ofile_.reset(nullptr);
    nWritten_ = 0;
    coordSetWriter::close(force);
}


void Foam::coordSetWriters::xmgraceWriter::beginTime(const Time& t)
{
    ofile_.reset(nullptr);
    nWritten_ = 0;
    coordSetWriter::beginTime(t);
}


void Foam::coordSetWriters::xmgraceWriter::beginTime(const instant& inst)
{
    ofile_.reset(nullptr);
    nWritten_ = 0;
    coordSetWriter::beginTime(inst);
}


void Foam::coordSetWriters::xmgraceWriter::endTime()
{
    ofile_.reset(nullptr);
    nWritten_ = 0;
    coordSetWriter::endTime();
}


template<class Type>
Foam::fileName Foam::coordSetWriters::xmgraceWriter::writeTemplate
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


    // Regular version

    const auto& coords = coords_[0];

    if (!ofile_)
    {
        // Field:
        // 1) rootdir/<TIME>/setName.agr
        // 2) rootdir/setName.agr

        const fileName outputFile = path();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        ofile_.reset(new OFstream(outputFile, streamOpt_));
        auto& os = ofile_();
        os.precision(precision_);

        // Preamble
        os  << "@g0 on" << nl
            << "@with g0" << nl
            << "@    title \"" << coords.name() << '"' << nl
            << "@    xaxis label \"" << coords.axis() << '"' << nl;

        nWritten_ = 0;  // Restarted
    }
    auto& os = ofile_();

    // Plot entry
    {
        os  << "@    s" << nWritten_
            << " legend \"" << fieldName << '"' << nl
            << "@target G0.S" << nWritten_ << nl;

        writeTable(os, coords, values, " \t");

        os  << '&' << nl;
        os  << "# end_data" << nl;
        ++nWritten_;
    }

    return ofile_().name();
}


template<class Type>
Foam::fileName Foam::coordSetWriters::xmgraceWriter::writeTemplate
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
defineCoordSetWriterWriteFields(Foam::coordSetWriters::xmgraceWriter);


// ************************************************************************* //
