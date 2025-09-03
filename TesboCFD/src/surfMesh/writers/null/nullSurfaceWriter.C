/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "nullSurfaceWriter.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(nullWriter);
    addToRunTimeSelectionTable(surfaceWriter, nullWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, nullWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::nullWriter::nullWriter()
:
    surfaceWriter()
{}


Foam::surfaceWriters::nullWriter::nullWriter(const dictionary& options)
:
    nullWriter()
{}


Foam::surfaceWriters::nullWriter::nullWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    nullWriter()
{}


Foam::surfaceWriters::nullWriter::nullWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    nullWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfaceWriters::nullWriter::needsUpdate() const
{
    return false;
}


bool Foam::surfaceWriters::nullWriter::wroteData() const
{
    return true;
}


bool Foam::surfaceWriters::nullWriter::enabled() const
{
    return false;
}


void Foam::surfaceWriters::nullWriter::setSurface
(
    const meshedSurf& surf,
    bool parallel
)
{}


void Foam::surfaceWriters::nullWriter::setSurface
(
    const pointField& points,
    const faceList& faces,
    bool parallel
)
{}


void Foam::surfaceWriters::nullWriter::open(const fileName& outputPath)
{}


Foam::fileName Foam::surfaceWriters::nullWriter::write()
{
    wroteGeom_ = true;
    return fileName::null;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::nullWriter);


// ************************************************************************* //
