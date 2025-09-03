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

#include "csvCoordSetWriter.H"
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
    defineTypeName(csvWriter);
    addToRunTimeSelectionTable(coordSetWriter, csvWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, csvWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "csvCoordSetWriterImpl.C"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::csvWriter::csvWriter()
:
    coordSetWriter(),
    streamOpt_(),
    precision_(IOstream::defaultPrecision())
{
    buffering_ = true;
}


Foam::coordSetWriters::csvWriter::csvWriter(const dictionary& options)
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


Foam::coordSetWriters::csvWriter::csvWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    csvWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::csvWriter::csvWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    csvWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::csvWriter::~csvWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSetWriters::csvWriter::buffering(const bool on)
{
    const bool old(buffering_);
    buffering_ = on;
    return old;
}


Foam::fileName Foam::coordSetWriters::csvWriter::path() const
{
    // Assume !useTracks_, otherwise too fragile

    // 1) rootdir/<TIME>/setName.csv
    // 2) rootdir/setName.csv

    return getExpectedPath("csv");
}


bool Foam::coordSetWriters::csvWriter::writeBuffered()
{
    if (coords_.empty())
    {
        clearBuffers();
        return false;
    }
    const auto& coords = coords_[0];


    DynamicList<word> headCols(3 + nDataColumns());

    do
    {
        if (coords.hasVectorAxis())
        {
            // x, y, z
            headCols.append("x");
            headCols.append("y");
            headCols.append("z");
        }
        else
        {
            headCols.append(coords.axis());
        }

        // label, scalar
        headCols.append(labelNames_);
        headCols.append(scalarNames_);

        // vector space
        #undef doLocalCode
        #define doLocalCode(Type)                                             \
                                                                              \
            for (const word& fldName : Type##Names_)                          \
            {                                                                 \
                for (direction d=0; d < pTraits<Type>::nComponents; ++d)      \
                {                                                             \
                    headCols.append(fldName + '_' + Foam::name(d));           \
                }                                                             \
            }

        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);
        #undef doLocalCode
    }
    while (false);


    // Field:
    // 1) rootdir/<TIME>/setName.csv
    // 2) rootdir/setName.csv

    fileName outputFile = path();

    if (!isDir(outputFile.path()))
    {
        mkDir(outputFile.path());
    }

    OFstream os(outputFile, streamOpt_);
    os.precision(precision_);

    writeLine(os, headCols, ",");

    writeBufferContents(os, coords, ",");

    clearBuffers();

    return true;
}


// * * * * * * * * * * * * * * * Implementation * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::coordSetWriters::csvWriter::writeTemplate
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
Foam::fileName Foam::coordSetWriters::csvWriter::writeTemplate
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
defineCoordSetWriterWriteFields(Foam::coordSetWriters::csvWriter);


// ************************************************************************* //
