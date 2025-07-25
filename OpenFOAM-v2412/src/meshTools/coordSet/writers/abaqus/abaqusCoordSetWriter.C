/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "abaqusCoordSetWriter.H"
#include "coordSet.H"
#include "IOmanip.H"
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
    defineTypeName(abaqusWriter);
    addToRunTimeSelectionTable(coordSetWriter, abaqusWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, abaqusWriter, wordDict);
}
}

const Foam::Enum<Foam::coordSetWriters::abaqusWriter::timeBase>
Foam::coordSetWriters::abaqusWriter::timeBaseNames_
({
    { timeBase::time, "time" },
    { timeBase::iter, "iteration" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
template<class Type>
static inline void putValue(Ostream& os, const Type& value, const int width)
{
    if (width) os << setw(width);
    os << value;
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::coordSetWriters::abaqusWriter::replaceUserEntries
(
    const string& str,
    const dictionary& vars
) const
{
    string result = str;

    const bool allowEnv = true;
    const bool allowEmpty = false;

    stringOps::inplaceExpand(result, vars, allowEnv, allowEmpty);

    return result;
}


void Foam::coordSetWriters::abaqusWriter::appendTimeName
(
    const word& fieldName,
    fileName& fName
) const
{
    if (useTimeDir())
    {
        return;
    }

    switch (timeBase_)
    {
        case timeBase::time:
        {
            fName.ext(timeName());
            break;
        }
        case timeBase::iter:
        {
            fName.ext(Foam::name(writeIndex_[fieldName]));
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << timeBaseNames_[timeBase_]
                << ". Available options: " << timeBaseNames_.sortedToc()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::abaqusWriter::abaqusWriter()
:
    coordSetWriter(),
    outputHeader_(),
    writeGeometry_(false),
    nullValue_(pTraits<scalar>::min),
    useLocalTimeDir_(coordSetWriter::useTimeDir()),
    timeBase_(timeBase::time),
    writeIndex_(0)
{}


Foam::coordSetWriters::abaqusWriter::abaqusWriter(const dictionary& options)
:
    coordSetWriter(options),
    outputHeader_(),
    writeGeometry_(false),
    nullValue_(pTraits<scalar>::min),
    useLocalTimeDir_(coordSetWriter::useTimeDir()),
    timeBase_(timeBase::time),
    writeIndex_(0)
{
    options.readIfPresent("header", outputHeader_);

    options.readIfPresent("useTimeDir", useLocalTimeDir_);

    if (!useLocalTimeDir_)
    {
        timeBaseNames_.readIfPresent("timeBase", options, timeBase_);
        options.readIfPresent("writeIndex", writeIndex_);
    }

    options.readIfPresent("writeGeometry", writeGeometry_);

    options.readIfPresent("nullValue", nullValue_);
}


Foam::coordSetWriters::abaqusWriter::abaqusWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    abaqusWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::abaqusWriter::abaqusWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    abaqusWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::abaqusWriter::~abaqusWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::coordSetWriters::abaqusWriter::path() const
{
    // 1) rootdir/<TIME>/setName.{inp}
    // 2) rootdir/setName.{inp}

    return getExpectedPath("inp");
}


void Foam::coordSetWriters::abaqusWriter::writeGeometry
(
    Ostream& os,
    label nTracks
) const
{
    if (!writeGeometry_ || coords_.empty())
    {
        return;
    }

    os  << "** Geometry" << nl
        << "**" << nl
        << "** Points" << nl
        << "**" << nl;

    // Write points
    label globalPointi = 1;
    for (const coordSet& coords : coords_)
    {
        for (const point& p : coords)
        {
            const point tp = p*geometryScale_;

            os  << globalPointi << ", "
                << tp[0] << ", " << tp[1] << ", " << tp[2]  << nl;

            ++globalPointi;
        }
    }

    if (nTracks)
    {
        WarningInFunction
            << "Tracks not implemented for " << typeName << endl;
    }

    wroteGeom_ = true;
}


template<class Type>
Foam::fileName Foam::coordSetWriters::abaqusWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& values
)
{
    // useTimeDir(options.getOrDefault("useTimeDir", true));
    useTimeDir(useLocalTimeDir_);

    // Note - invalid samples are set to pTraits<Type>::max in sampledSetsImpl.C

    checkOpen();
    if (coords_.empty())
    {
        return fileName::null;
    }

    fileName outputFile = path();

    if (!wroteGeom_)
    {
        if (!writeIndex_.insert(fieldName, 0))
        {
            ++writeIndex_[fieldName];
        }

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        if (writeGeometry_)
        {
            const word geomName("geometry");
            if (!writeIndex_.insert(geomName, 0))
            {
                ++writeIndex_[geomName];
            }

            fileName geomFileName(outputFile.lessExt() + ".inp");

            appendTimeName("geometry", geomFileName);

            if (verbose_)
            {
                Info<< "Writing abaqus geometry to " << geomFileName << endl;
            }

            OFstream osGeom(geomFileName);

            writeGeometry(osGeom, (useTracks_ ? coords_.size() : 0));
        }

        fileName fieldFileName(outputFile.lessExt() + ".inp_" + fieldName);
        appendTimeName(fieldName, fieldFileName);

        if (verbose_)
        {
            Info<< "Writing field data to " << fieldFileName << endl;
        }

        OFstream os(fieldFileName);

        if (!outputHeader_.empty())
        {
            dictionary vars;
            vars.add("TIME", timeName());
            vars.add("FIELD_NAME", fieldName);
            vars.add("FILE_NAME", fieldFileName);

            for (const auto& s : outputHeader_)
            {
                os  << replaceUserEntries(s, vars).c_str() << nl;
            }
        }
        else
        {
            os  << "** OpenFOAM " << fieldFileName.nameLessExt() << nl
                << "** Project " << outputFile << nl
                << "** Field=" << fieldName << " Time=" << timeName() << nl;
        }

        tmp<Field<Type>> tfield(values);
        tfield = adjustFieldTemplate(fieldName, tfield);
        const auto& field = tfield();

        forAll(field, samplei)
        {
            os  << (samplei+1);
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
            {
                // Work-around to set null values - set to pTraits<Type>::max
                // by default
                scalar s = component(field[samplei], cmpt);
                if (s > 0.5*pTraits<scalar>::max)
                {
                    s = nullValue_;
                }

                os  << ", " << s;
            }
            os  << nl;
        }
    }

    return outputFile;
}


template<class Type>
Foam::fileName Foam::coordSetWriters::abaqusWriter::writeTemplate
(
    const word& fieldName,
    const List<Field<Type>>& fieldValues
)
{
    // Track writing

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
            Info<< "Writing abaqus geometry to " << outputFile << endl;
        }

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream osGeom(outputFile.lessExt() + "." + fieldName + ".inp");

        osGeom
            << "** Geometry" << nl
            << "**" << nl
            << "** Points" << nl
            << "**" << nl;

        writeGeometry(osGeom, coords_.size());
    }

    return outputFile;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineCoordSetWriterWriteFields(Foam::coordSetWriters::abaqusWriter);


// ************************************************************************* //
