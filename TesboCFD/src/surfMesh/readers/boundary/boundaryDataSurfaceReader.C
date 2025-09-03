/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "boundaryDataSurfaceReader.H"
#include "rawIOField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundaryDataSurfaceReader, 0);
    addToRunTimeSelectionTable
    (
        surfaceReader,
        boundaryDataSurfaceReader,
        fileName
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::pointField Foam::boundaryDataSurfaceReader::readPoints
(
    const Time& runTime,
    const fileName& baseDir,
    const word& pointsName
)
{
    fileName pointsFile
    (
        baseDir / (pointsName.empty() ? "points" : pointsName)
    );
    pointsFile.toAbsolute();

    IOobject io
    (
        pointsFile,             // absolute path
        runTime,
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER,
        true                    // global object (currently not used)
    );

    DebugInfo
        << "File: " << io.objectPath() << endl;

    // Read data (no average value!)
    rawIOField<point> rawData(io);

    pointField points(std::move(rawData.field()));

    DebugInfo
        << "File: " << io.objectPath()
        << " " << points.size() << " points" << endl;

    return points;
}


Foam::pointField Foam::boundaryDataSurfaceReader::readPoints
(
    const fileName& dirName,
    const word& pointsName
)
{
    refPtr<Time> timePtr(Time::NewGlobalTime());

    return readPoints(*timePtr, dirName, pointsName);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundaryDataSurfaceReader::readCase()
{
    DebugInFunction << endl;

    timeValues_ = TimePaths::findTimes(baseDir_);
}


void Foam::boundaryDataSurfaceReader::readGeometry
(
    meshedSurface& surf,
    const label timeIndex
)
{
    surf.clear();

    pointField points(this->readPoints(baseDir_, pointsName_));

    surf = meshedSurface(std::move(points), faceList());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::boundaryDataSurfaceReader::readFieldTemplate
(
    const label timeIndex,
    const label fieldIndex
) const
{
    Type dummyAvg;

    return readField<Type>
    (
        baseDir_,
        timeValues_[timeIndex],
        fieldNames_[fieldIndex],
        dummyAvg
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryDataSurfaceReader::boundaryDataSurfaceReader
(
    const fileName& fName,
    const word& pointsName
)
:
    boundaryDataSurfaceReader(fName, dictionary(), pointsName)
{}


Foam::boundaryDataSurfaceReader::boundaryDataSurfaceReader
(
    const fileName& fName,
    const dictionary& options,
    const word& pointsName
)
:
    surfaceReader(fName, options),
    baseDir_(fName.path()),
    pointsName_(pointsName),
    timeValues_(),
    fieldNames_(),
    surfPtr_(nullptr)
{
    options.readIfPresent("points", pointsName_);

    baseDir_.toAbsolute();
    debug = 1;
    DebugInFunction << endl;
    Info<< "create with " << baseDir_ << endl;

    readCase();
}


// * * * * * * * * * * * * * Public Member Functions   * * * * * * * * * * * //

const Foam::meshedSurface& Foam::boundaryDataSurfaceReader::geometry
(
    const label timeIndex
)
{
    DebugInFunction << endl;

    if (!surfPtr_)
    {
        surfPtr_.reset(new meshedSurface);
        readGeometry(*surfPtr_, timeIndex);
    }

    return *surfPtr_;
}


Foam::instantList Foam::boundaryDataSurfaceReader::times() const
{
    return timeValues_;
}


Foam::wordList Foam::boundaryDataSurfaceReader::fieldNames
(
    const label timeIndex
) const
{
    if (timeValues_.empty() || timeIndex >= timeValues_.size())
    {
        fieldNames_.clear();
        return wordList();
    }

    fileNameList items =
        fileHandler().readDir(baseDir_/timeValues_[timeIndex].name());

    fieldNames_.resize_nocopy(items.size());

    std::transform
    (
        items.begin(),
        items.end(),
        fieldNames_.begin(),
        [](const fileName& f) { return word(f); }
    );

    Foam::sort(fieldNames_);

    return fieldNames_;
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::boundaryDataSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const scalar& refValue
) const
{
    return readFieldTemplate<scalar>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::vector>> Foam::boundaryDataSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const vector& refValue
) const
{
    return readFieldTemplate<vector>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::boundaryDataSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const sphericalTensor& refValue
) const
{
    return readFieldTemplate<sphericalTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::symmTensor>> Foam::boundaryDataSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const symmTensor& refValue
) const
{
    return readFieldTemplate<symmTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::tensor>> Foam::boundaryDataSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const tensor& refValue
) const
{
    return readFieldTemplate<tensor>(timeIndex, fieldIndex);
}


// ************************************************************************* //
