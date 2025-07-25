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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::boundaryDataSurfaceReader::readField
(
    const Time& runTime,
    const fileName& baseDir,
    const instant& timeDir,
    const word& fieldName,
    Type& avg
)
{
    // Reread values and interpolate
    fileName valuesFile(baseDir/timeDir.name()/fieldName);
    valuesFile.toAbsolute();

    IOobject io
    (
        valuesFile,   // absolute path
        runTime,
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER,
        true                // global object (currently not used)
    );

    DebugInfo
        << "File: " << io.objectPath() << endl;

    // Read data (TDB: setAverage)
    rawIOField<Type> rawData(io, IOobjectOption::LAZY_READ);

    if (rawData.hasAverage())
    {
        avg = rawData.average();
    }

    DebugInfo
        << "File: " << io.objectPath()
        << " " << rawData.size() << " values" << endl;

    return tmp<Field<Type>>::New(std::move(rawData.field()));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::boundaryDataSurfaceReader::readField
(
    const fileName& baseDir,
    const instant& timeDir,
    const word& fieldName,
    Type& avg
)
{
    refPtr<Time> timePtr(Time::NewGlobalTime());

    return readField<Type>(*timePtr, baseDir, timeDir, fieldName, avg);
}


// ************************************************************************* //
