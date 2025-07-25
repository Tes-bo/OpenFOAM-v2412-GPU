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

#include "nullCoordSetWriter.H"
#include "coordSetWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{
    defineTypeName(nullWriter);
    addToRunTimeSelectionTable(coordSetWriter, nullWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, nullWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::nullWriter::nullWriter()
:
    coordSetWriter()
{}


Foam::coordSetWriters::nullWriter::nullWriter(const dictionary& options)
:
    nullWriter()
{}


Foam::coordSetWriters::nullWriter::nullWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    nullWriter()
{}


Foam::coordSetWriters::nullWriter::nullWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    nullWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::nullWriter::~nullWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSetWriters::nullWriter::enabled() const
{
    return false;
}


bool Foam::coordSetWriters::nullWriter::buffering() const
{
    return false;
}


bool Foam::coordSetWriters::nullWriter::needsUpdate() const
{
    return false;
}


bool Foam::coordSetWriters::nullWriter::wroteData() const
{
    return true;
}


void Foam::coordSetWriters::nullWriter::setCoordinates
(
    const coordSet* coords
)
{}


void Foam::coordSetWriters::nullWriter::setCoordinates
(
    const coordSet& coords
)
{}


void Foam::coordSetWriters::nullWriter::setTracks
(
    const UPtrList<coordSet>& tracks
)
{}


Foam::fileName Foam::coordSetWriters::nullWriter::path() const
{
    return fileName();
}


void Foam::coordSetWriters::nullWriter::open(const fileName& outputPath)
{}


// Foam::fileName Foam::coordSetWriters::nullWriter::write()
// {
//     wroteGeom_ = true;
//     return fileName::null;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineCoordSetWriterWriteFields(Foam::coordSetWriters::nullWriter);


// ************************************************************************* //
