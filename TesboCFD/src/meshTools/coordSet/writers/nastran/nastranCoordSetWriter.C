/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "nastranCoordSetWriter.H"
#include "coordSet.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "coordSetWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{
    defineTypeName(nastranWriter);
    addToRunTimeSelectionTable(coordSetWriter, nastranWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, nastranWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static inline void putValue(Ostream& os, const Type& value, const int width)
{
    if (width) os << setw(width);
    os << value;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Ostream& Foam::coordSetWriters::nastranWriter::writeKeyword
(
    Ostream& os,
    const word& keyword
) const
{
    return fileFormats::NASCore::writeKeyword(os, keyword, writeFormat_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::nastranWriter::nastranWriter()
:
    coordSetWriter(),
    writeFormat_(fieldFormat::FREE),
    separator_()
{
    if (writeFormat_ == fieldFormat::FREE)
    {
        separator_ = ",";
    }
}


Foam::coordSetWriters::nastranWriter::nastranWriter(const dictionary& options)
:
    coordSetWriter(options),
    writeFormat_
    (
        fileFormats::NASCore::fieldFormatNames.getOrDefault
        (
            "format",
            options,
            fieldFormat::FREE
        )
    ),
    separator_()
{
    if (writeFormat_ == fieldFormat::FREE)
    {
        separator_ = ",";
    }
}


Foam::coordSetWriters::nastranWriter::nastranWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    nastranWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::nastranWriter::nastranWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    nastranWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::nastranWriter::~nastranWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::coordSetWriters::nastranWriter::path() const
{
    // 1) rootdir/<TIME>/setName.{nas}
    // 2) rootdir/setName.{nas}

    return getExpectedPath("nas");
}


void Foam::coordSetWriters::nastranWriter::writeGeometry
(
    Ostream& os,
    label nTracks
) const
{
    if (coords_.empty())
    {
        return;
    }

    // Field width (SHORT, LONG formats)
    const int width =
    (
        writeFormat_ == fieldFormat::SHORT ? 8
      : writeFormat_ == fieldFormat::LONG ? 16
      : 0
    );

    // Separator char (FREE format)
    const char sep = (writeFormat_ == fieldFormat::FREE ? ',' : '\0');

    // Write points
    os  << '$' << nl
        << "$ Points" << nl
        << '$' << nl;

    label globalPointi = 0;
    for (const coordSet& coords : coords_)
    {
        for (const point& p : coords)
        {
            fileFormats::NASCore::writeCoord
            (
                os, p, globalPointi, writeFormat_
            );
            ++globalPointi;
        }
    }

    if (nTracks)
    {
        // Write ids of track points to file
        globalPointi = 0;
        label globalEdgei = 0;

        for (label tracki = 0; tracki < nTracks; ++tracki)
        {
            const label nEdges = (coords_[tracki].size() - 1);

            for (label edgei = 0; edgei < nEdges; ++edgei)
            {
                writeKeyword(os, "PLOTEL");
                if (sep) os << sep;

                putValue(os, globalEdgei+1, width);  // Edge id
                if (sep) os << sep;

                putValue(os, globalPointi+1, width);
                if (sep) os << sep;

                putValue(os, globalPointi+2, width);
                os << nl;

                ++globalEdgei;
                ++globalPointi;
            }
        }
    }

    wroteGeom_ = true;
}


template<class Type>
Foam::fileName Foam::coordSetWriters::nastranWriter::writeTemplate
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

    fileName outputFile = path();

    if (!wroteGeom_)
    {
        if (verbose_)
        {
            Info<< "Writing nastran geometry to " << outputFile << endl;
        }

        if (!Foam::isDir(outputFile.path()))
        {
            Foam::mkDir(outputFile.path());
        }

        OFstream os(IOstreamOption::ATOMIC, outputFile);
        fileFormats::NASCore::setPrecision(os, writeFormat_);

        os  << "TITLE=OpenFOAM " << outputFile.stem()
            << " geometry" << nl
            << "BEGIN BULK" << nl;

        writeGeometry(os, (useTracks_ ? coords_.size() : 0));

        os << "ENDDATA" << nl;
    }

    return outputFile;
}


template<class Type>
Foam::fileName Foam::coordSetWriters::nastranWriter::writeTemplate
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

    fileName outputFile = path();

    if (!wroteGeom_)
    {
        if (verbose_)
        {
            Info<< "Writing nastran geometry to " << outputFile << endl;
        }

        if (!Foam::isDir(outputFile.path()))
        {
            Foam::mkDir(outputFile.path());
        }

        OFstream os(IOstreamOption::ATOMIC, outputFile);
        fileFormats::NASCore::setPrecision(os, writeFormat_);

        os  << "TITLE=OpenFOAM " << outputFile.stem()
            << " geometry" << nl
            << "BEGIN BULK" << nl;

        writeGeometry(os, coords_.size());

        os << "ENDDATA" << nl;
    }

    return outputFile;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineCoordSetWriterWriteFields(Foam::coordSetWriters::nastranWriter);


// ************************************************************************* //
