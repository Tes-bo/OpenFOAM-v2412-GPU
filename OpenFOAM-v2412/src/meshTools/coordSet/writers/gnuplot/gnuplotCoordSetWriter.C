/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "gnuplotCoordSetWriter.H"
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
    defineTypeName(gnuplotWriter);
    addToRunTimeSelectionTable(coordSetWriter, gnuplotWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, gnuplotWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "gnuplotCoordSetWriterImpl.C"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::gnuplotWriter::gnuplotWriter()
:
    coordSetWriter(),
    streamOpt_(),
    precision_(IOstream::defaultPrecision())
{
    buffering_ = true;
}


Foam::coordSetWriters::gnuplotWriter::gnuplotWriter(const dictionary& options)
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


Foam::coordSetWriters::gnuplotWriter::gnuplotWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    gnuplotWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::gnuplotWriter::gnuplotWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    gnuplotWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::gnuplotWriter::~gnuplotWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSetWriters::gnuplotWriter::buffering(const bool on)
{
    const bool old(buffering_);
    buffering_ = on;
    return old;
}


Foam::fileName Foam::coordSetWriters::gnuplotWriter::path() const
{
    // 1) rootdir/<TIME>/setName.{gplt}
    // 2) rootdir/setName.{gplt}

    return getExpectedPath("gplt");
}


bool Foam::coordSetWriters::gnuplotWriter::writeBuffered()
{
    if (coords_.empty())
    {
        clearBuffers();
        return false;
    }

    // Field:
    // 1) rootdir/<TIME>/setName.gplt
    // 2) rootdir/setName.gplt

    fileName outputFile = path();

    if (!isDir(outputFile.path()))
    {
        mkDir(outputFile.path());
    }

    OFstream os(outputFile, streamOpt_);
    os.precision(precision_);

    os  << "set term pngcairo" << nl
        << "set output \"" << outputFile.stem() << ".png\"" << nl;

    label nplots = 0;
    do
    {
        #undef doLocalCode
        #define doLocalCode(Type)                                             \
        for (const word& fldName : Type##Names_)                              \
        {                                                                     \
            os  << (nplots++ ? ", \\" : "plot \\") << nl;                     \
            os  << "  '-' title \"" << fldName << "\" with lines";            \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);
        #undef doLocalCode
    }
    while (false);

    os  << nl << nl;

    if (nplots)
    {
        #undef doLocalCode
        #define doLocalCode(Type)                                             \
        for (const Field<Type>& fld : Type##Fields_)                          \
        {                                                                     \
            writeTable(os, coords_[0], fld, " \t");                           \
            os  << "end_data" << nl << nl;                                    \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);
        #undef doLocalCode
    }

    os  << "# end plot" << nl;

    clearBuffers();

    return true;
}


template<class Type>
Foam::fileName Foam::coordSetWriters::gnuplotWriter::writeTemplate
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
Foam::fileName Foam::coordSetWriters::gnuplotWriter::writeTemplate
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

    UPtrList<const Field<Type>> fieldPtrs(repackageFields(fieldValues));
    return writeTemplate(fieldName, fieldPtrs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineCoordSetWriterWriteFields(Foam::coordSetWriters::gnuplotWriter);


// ************************************************************************* //
