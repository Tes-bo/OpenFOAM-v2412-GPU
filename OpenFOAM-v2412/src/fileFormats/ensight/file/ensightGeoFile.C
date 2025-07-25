/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "ensightGeoFile.H"
#include "foamVersion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightGeoFile::ensightGeoFile
(
    IOstreamOption::appendType append,
    const fileName& pathname,
    IOstreamOption::streamFormat fmt
)
:
    ensightFile(append, pathname, fmt)
{
    if (!OFstream::is_appending())
    {
        writeBinaryHeader();  // Mandatory for geometry files
    }
}


Foam::ensightGeoFile::ensightGeoFile
(
    IOstreamOption::appendType append,
    const fileName& path,
    const fileName& name,
    IOstreamOption::streamFormat fmt
)
:
    ensightFile(append, path, name, fmt)
{
    if (!OFstream::is_appending())
    {
        writeBinaryHeader();  // Mandatory for geometry files
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::ensightGeoFile::writeKeyword(const keyType& key)
{
    writeString(key);
    newline();

    return *this;
}


//
// Convenience Output Methods
//

void Foam::ensightGeoFile::beginGeometry()
{
    // Description line 1
    writeString("Ensight Geometry File");
    newline();

    // Description line 2
    writeString("Written by OpenFOAM " + std::to_string(foamVersion::api));
    newline();

    writeString("node id assign");
    newline();

    writeString("element id assign");
    newline();
}


// ************************************************************************* //
