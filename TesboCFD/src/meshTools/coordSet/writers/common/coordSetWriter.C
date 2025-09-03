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

#include "coordSet.H"
#include "coordSetWriter.H"
#include "fileFormats.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordSetWriter, 0);
    defineRunTimeSelectionTable(coordSetWriter, word);
    defineRunTimeSelectionTable(coordSetWriter, wordDict);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::coordSetWriter::suffix
(
    const word& fldName,
    const word& fileExt
)
{
    word result;

    if (!fldName.empty())
    {
        result += '_' + fldName;
    }

    return result.ext(fileExt);
}


Foam::word Foam::coordSetWriter::suffix
(
    const wordList& fieldNames,
    const word& fileExt
)
{
    word result;

    for (const word& fldName : fieldNames)
    {
        if (!fldName.empty())
        {
            result += '_' + fldName;
        }
    }

    return result.ext(fileExt);
}


Foam::dictionary Foam::coordSetWriter::formatOptions
(
    const dictionary& dict,
    const word& formatName,
    const word& entryName
)
{
    return fileFormats::getFormatOptions(dict, formatName, entryName);
}


Foam::dictionary Foam::coordSetWriter::formatOptions
(
    const dictionary& dict,
    const dictionary& setDict,
    const word& formatName,
    const word& entryName
)
{
    return fileFormats::getFormatOptions(dict, setDict, formatName, entryName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriter::coordSetWriter()
:
    coords_(),
    trackTimes_(),
    upToDate_(false),
    wroteGeom_(false),
/// parallel_(true),
    buffering_(false),
    useTracks_(false),
    useTimeDir_(false),
    verbose_(false),
    nFields_(0),
    currTime_(),
    outputPath_(),
    geometryScale_(1),
    geometryCentre_(Zero),
    geometryTransform_(),
    fieldLevel_(),
    fieldScale_()
{}


Foam::coordSetWriter::coordSetWriter(const dictionary& options)
:
    coordSetWriter()
{
    options.readIfPresent("verbose", verbose_);

    geometryScale_ = 1;
    geometryCentre_ = Zero;
    geometryTransform_.clear();

    options.readIfPresent("scale", geometryScale_);

    // Optional cartesian coordinate system transform
    const auto* dictptr = options.findDict("transform", keyType::LITERAL);

    if (dictptr)
    {
        dictptr->readIfPresent("rotationCentre", geometryCentre_);

        // 'origin' is optional within sub-dictionary
        geometryTransform_ =
            coordSystem::cartesian(*dictptr, IOobjectOption::LAZY_READ);
    }

    fieldLevel_ = options.subOrEmptyDict("fieldLevel");
    fieldScale_ = options.subOrEmptyDict("fieldScale");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriter::~coordSetWriter()
{
    close();
}


// * * * * * * * * * * * * * * * * * Controls  * * * * * * * * * * * * * * * //

void Foam::coordSetWriter::setTime(const instant& inst)
{
    currTime_ = inst;
}


void Foam::coordSetWriter::setTime(scalar timeValue)
{
    currTime_ = instant(timeValue);
}


void Foam::coordSetWriter::setTime(scalar timeValue, const word& timeName)
{
    currTime_.value() = timeValue;
    currTime_.name() = timeName;
}


void Foam::coordSetWriter::unsetTime()
{
    currTime_.value() = 0;
    currTime_.name().clear();
}


void Foam::coordSetWriter::beginTime(const Time& t)
{
    setTime(t.value(), t.timeName());
}


void Foam::coordSetWriter::beginTime(const instant& inst)
{
    setTime(inst);
}


void Foam::coordSetWriter::endTime()
{
    // Flush bufferred data
    if (nDataColumns())
    {
        writeBuffered();
    }
    clearBuffers();

    unsetTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coordSetWriter::open(const fileName& outputPath)
{
    outputPath_ = outputPath;
    wroteGeom_ = false;
}


void Foam::coordSetWriter::open
(
    const coordSet& coords,
    const fileName& outputPath
)
{
    close();
    setCoordinates(coords);
    open(outputPath);
}


void Foam::coordSetWriter::open
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath
)
{
    close();
    setTracks(tracks);
    open(outputPath);
}


void Foam::coordSetWriter::close(bool force)
{
    if (nDataColumns())
    {
        if (verbose_) Info<< "Flush buffered data:" << nl;
        writeBuffered();
    }
    clearBuffers();

    outputPath_.clear();
    wroteGeom_ = false;

    if (force)
    {
        coords_.clear();
        trackTimes_.clear();
    }
}


void Foam::coordSetWriter::clear()
{
    close();
    expire();
    coords_.clear();
    trackTimes_.clear();
    clearBuffers();  // Reset any buffering
}


void Foam::coordSetWriter::setCoordinates(const coordSet* coords)
{
    expire();
    clearBuffers();  // Reset any buffering

    if (coords)
    {
        coords_.resize(1);
        coords_.set(0, coords);
    }
    else
    {
        coords_.clear();
    }
    trackTimes_.clear();
}


void Foam::coordSetWriter::setCoordinates(const coordSet& coords)
{
    setCoordinates(&coords);
}


void Foam::coordSetWriter::setTracks(const UPtrList<coordSet>& tracks)
{
    expire();
    clearBuffers();  // Reset any buffering

    // Shallow copy (pointers)

    coords_.resize(tracks.size());
    forAll(coords_, tracki)
    {
        coords_.set(tracki, tracks.get(tracki));
    }
    trackTimes_.clear();
    useTracks_ = true;
}


void Foam::coordSetWriter::setTrackTimes(const UList<scalarField>& times)
{
    if (times.size() == coords_.size())
    {
        trackTimes_ = times;
    }
    else
    {
        trackTimes_.clear();
    }
}


Foam::label Foam::coordSetWriter::numPoints() const
{
    label nTotal = 0;

    forAll(coords_, tracki)
    {
        const auto* ptr = coords_.get(tracki);
        if (ptr) nTotal += ptr->size();
    }

    return nTotal;
}


Foam::label Foam::coordSetWriter::numTracks() const
{
    return coords_.size();
}


bool Foam::coordSetWriter::needsUpdate() const
{
    return !upToDate_;
}


bool Foam::coordSetWriter::wroteData() const
{
    return wroteGeom_;
}


bool Foam::coordSetWriter::expire()
{
    const bool changed = upToDate_;

    upToDate_ = false;
    wroteGeom_ = false;
    coords_.clear();
    trackTimes_.clear();

    // Field count (nFields_) is a different type of accounting
    // and is unaffected by geometry changes

    return changed;
}


bool Foam::coordSetWriter::hasCoords() const
{
    return !coords_.empty();
}


bool Foam::coordSetWriter::empty() const
{
    return coords_.empty();
}

Foam::fileName Foam::coordSetWriter::getExpectedPath
(
    const word& fileExt
) const
{
    fileName file;

    if (!outputPath_.empty())
    {
        if (useTimeDir() && !timeName().empty())
        {
            // Splice in time-directory
            file = outputPath_.path() / timeName() / outputPath_.name();
        }
        else
        {
            file = outputPath_;
        }

        file.ext(fileExt);   // Append extension - can also be empty
    }

    return file;
}


Foam::fileName Foam::coordSetWriter::getFieldPrefixedPath
(
    const word& fieldName,
    const word& fileExt
) const
{
    if (outputPath_.empty() || fieldName.empty())
    {
        return getExpectedPath(fileExt);
    }

    // Field:  rootdir/<TIME>/<field>_NAME.ext

    fileName file;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        file = outputPath_.path() / timeName();
    }
    else
    {
        file = outputPath_.path();
    }

    // Append <field>_NAME.EXT
    file /= (fieldName + '_' + outputPath_.name());
    file.ext(fileExt);   // Append extension - can also be empty

    return file;
}


void Foam::coordSetWriter::checkOpen() const
{
    if (!is_open())
    {
        FatalErrorInFunction
            << type() << " : Attempted to write without a path" << nl
            << exit(FatalError);
    }
}


bool Foam::coordSetWriter::merge() const
{
    bool changed = false;

    // Possible future requirement...
    //
    // if (parallel_ && Pstream::parRun() && !upToDate_)
    // {
    //     changed = merged_.merge(coords_, mergeDim_);
    // }
    upToDate_ = true;

    if (changed)
    {
        wroteGeom_ = false;
    }

    return changed;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<coordSetWriter>& iproxy
)
{
    const auto& w = *iproxy;

    os  << "coordSetWriter:"
        << " upToDate: " << w.upToDate_
        << " nFields: " << w.nFields_
        << " time: " << w.currTime_
        << " path: " << w.outputPath_ << endl;

    return os;
}


// ************************************************************************* //
