/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "surfaceReader.H"
#include "fileFormats.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceReader, 0);
    defineRunTimeSelectionTable(surfaceReader, fileName);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::dictionary Foam::surfaceReader::formatOptions
(
    const dictionary& dict,
    const word& formatName,
    const word& entryName
)
{
    return fileFormats::getFormatOptions(dict, formatName, entryName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceReader::surfaceReader
(
    const fileName& fName
)
:
    fileName_(fName)
{}


Foam::surfaceReader::surfaceReader
(
    const fileName& fName,
    const dictionary& options
)
:
    surfaceReader(fName)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceReader>
Foam::surfaceReader::New
(
    const word& readerType,
    const fileName& fName,
    const dictionary& options
)
{
    auto* ctorPtr = fileNameConstructorTable(readerType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "reader",
            readerType,
            *fileNameConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<surfaceReader>(ctorPtr(fName, options));
}


// ************************************************************************* //
