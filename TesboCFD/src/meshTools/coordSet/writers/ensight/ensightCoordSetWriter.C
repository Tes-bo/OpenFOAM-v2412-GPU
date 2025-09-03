/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "ensightCoordSetWriter.H"
#include "coordSet.H"
#include "IOmanip.H"
#include "ensightCase.H"
#include "ensightGeoFile.H"
#include "ensightOutput.H"
#include "ensightPTraits.H"
#include "coordSetWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{
    defineTypeName(ensightWriter);
    addToRunTimeSelectionTable(coordSetWriter, ensightWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, ensightWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static void writeTrackField
(
    ensightFile& os,
    const UPtrList<const Field<Type>>& fieldPtrs
)
{
    // Write field (serial only)
    os.write(ensightPTraits<Type>::typeName);
    os.newline();

    forAll(fieldPtrs, tracki)
    {
        const Field<Type>& fld = fieldPtrs[tracki];

        os.beginPart(tracki);  // Part index (0-based)

        // Skip if empty. This check is probably redundant
        if (fld.empty())
        {
            continue;
        }

        // Write as coordinates data
        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightFile::coordinates,
            fld,
            false  /* serial only! */
        );
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordSetWriters::ensightWriter::writeGeometry
(
    ensightGeoFile& os,
    elemOutputType elemOutput
) const
{
    // Writing tracks as x/y/z coordinates, optionally with points
    // or points/lines as elements.
    //
    // The requirements are so basic that they do not warrant an
    // ensightPart treatment at all.

    forAll(coords_, tracki)
    {
        const auto& coords = coords_[tracki];
        const label nPoints = coords.size();

        word partName("track" + Foam::name(tracki));
        if (coords_.size() == 1 && elemOutputType::WRITE_LINES != elemOutput)
        {
            partName = "sampled";
        }

        ensightOutput::Detail::writeCoordinates
        (
            os,
            tracki,   // Part index (0-based)
            partName,
            nPoints,
            static_cast<const pointField&>(coords),
            false     /* serial only! */
        );

        if (elemOutputType::WRITE_POINTS == elemOutput)
        {
            if (nPoints)
            {
                os.writeKeyword("point");
                os.write(nPoints);
                os.newline();
                for (label pointi = 0; pointi < nPoints; ++pointi)
                {
                    os.write(pointi+1);  // From 0-based to 1-based index
                    os.newline();
                }
            }
        }
        if (elemOutputType::WRITE_LINES == elemOutput)
        {
            const label nLines = (nPoints-1);
            if (nPoints == 1)
            {
                os.writeKeyword("point");
                os.write(nPoints);
                os.newline();
                for (label pointi = 0; pointi < nPoints; ++pointi)
                {
                    os.write(pointi+1);  // From 0-based to 1-based index
                    os.newline();
                }
            }
            else if (nLines > 0)
            {
                os.writeKeyword("bar2");
                os.write(nLines);
                os.newline();
                for (label pointi = 0; pointi < nLines; ++pointi)
                {
                    os.write(pointi+1);  // From 0-based to 1-based index
                    os.write(pointi+2);
                    os.newline();
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::ensightWriter::ensightWriter()
:
    coordSetWriter(),
    caseOpts_(IOstreamOption::BINARY),
    collateTimes_(true),
    caching_("fieldsDict")  // Historic name
{}


Foam::coordSetWriters::ensightWriter::ensightWriter(const dictionary& options)
:
    coordSetWriter(options),
    caseOpts_("format", options, IOstreamOption::BINARY),
    collateTimes_(options.getOrDefault("collateTimes", true)),
    caching_("fieldsDict")  // Historic name
{
    caseOpts_.timeFormat("timeFormat", options);
    caseOpts_.timePrecision("timePrecision", options);
}


Foam::coordSetWriters::ensightWriter::ensightWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    ensightWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::ensightWriter::ensightWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    ensightWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::ensightWriter::~ensightWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::coordSetWriters::ensightWriter::path() const
{
    // Assume collateTimes == true, otherwise too fragile

    // Collated
    // ========
    // CaseFile:  rootdir/NAME/NAME.case
    // Geometry:  rootdir/NAME/data/<index>/geometry
    // Field:     rootdir/NAME/data/<index>/field

    if (!outputPath_.empty())
    {
        return outputPath_ / (ensight::FileName(outputPath_.name()) + ".case");
    }

    return fileName();
}


void Foam::coordSetWriters::ensightWriter::close(const bool force)
{
    caching_.clear();
    coordSetWriter::close(force);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Field writing implementations

#include "ensightCoordSetWriterCollated.C"
#include "ensightCoordSetWriterUncollated.C"

template<class Type>
Foam::fileName Foam::coordSetWriters::ensightWriter::writeTemplate
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
    if (coords_.size() != 1)
    {
        FatalErrorInFunction
            << "Attempted to write field: " << fieldName
            << " (" << 1 << " entries) for "
            << coords_.size() << " sets" << nl
            << exit(FatalError);
    }

    UPtrList<const Field<Type>> fieldPtrs(repackageFields(values));

    elemOutputType elemOutput =
    (
        useTracks_
      ? elemOutputType::WRITE_LINES
      : elemOutputType::NO_ELEMENTS
    );

    if (collateTimes_)
    {
        return writeCollated(fieldName, fieldPtrs, elemOutput);
    }
    else
    {
        return writeUncollated(fieldName, fieldPtrs, elemOutput);
    }
}


template<class Type>
Foam::fileName Foam::coordSetWriters::ensightWriter::writeTemplate
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
    if (coords_.size() != fieldValues.size())
    {
        FatalErrorInFunction
            << "Attempted to write field: " << fieldName
            << " (" << fieldValues.size() << " entries) for "
            << coords_.size() << " sets" << nl
            << exit(FatalError);
    }

    UPtrList<const Field<Type>> fieldPtrs(repackageFields(fieldValues));

    if (collateTimes_)
    {
        return writeCollated
        (
            fieldName,
            fieldPtrs,
            elemOutputType::WRITE_LINES
        );
    }
    else
    {
        return writeUncollated
        (
            fieldName,
            fieldPtrs,
            elemOutputType::WRITE_LINES
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineCoordSetWriterWriteFields(Foam::coordSetWriters::ensightWriter);


// ************************************************************************* //
