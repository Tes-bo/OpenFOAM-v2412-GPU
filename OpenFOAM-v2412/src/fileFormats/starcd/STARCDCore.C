/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "STARCDCore.H"
#include "DynamicList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fileFormats::STARCDCore::fileHeader
>
Foam::fileFormats::STARCDCore::fileHeaders_
({
    { fileHeader::HEADER_CEL, "PROSTAR_CELL" },
    { fileHeader::HEADER_VRT, "PROSTAR_VERTEX" },
    { fileHeader::HEADER_BND, "PROSTAR_BOUNDARY" },
});


const Foam::Enum
<
    Foam::fileFormats::STARCDCore::fileExt
>
Foam::fileFormats::STARCDCore::fileExtensions_
({
    { fileExt::CEL_FILE, "cel" },
    { fileExt::VRT_FILE, "vrt" },
    { fileExt::BND_FILE, "bnd" },
    { fileExt::INP_FILE, "inp" },
});


const char* const Foam::fileFormats::STARCDCore::defaultBoundaryName =
    "Default_Boundary_Region";

const char* const Foam::fileFormats::STARCDCore::defaultSolidBoundaryName =
    "Default_Boundary_Solid";


const Foam::Map<Foam::FixedList<int, 6>>
Foam::fileFormats::STARCDCore::foamToStarFaceAddr =
{
    { starcdHex,   { 4, 5, 2, 3, 0, 1 } },
    { starcdPrism, { 0, 1, 4, 5, 2, -1 } },
    { starcdTet,   { 5, 4, 2, 0, -1, -1 } },
    { starcdPyr,   { 0, 4, 3, 5, 2, -1 } }
};


const Foam::Map<Foam::FixedList<int, 6>>
Foam::fileFormats::STARCDCore::starToFoamFaceAddr =
{
    { starcdHex,   { 4, 5, 2, 3, 0, 1 } },
    { starcdPrism, { 0, 1, 4, -1, 2, 3 } },
    { starcdTet,   { 3, -1, 2, -1, 1, 0 } },
    { starcdPyr,   { 0, -1, 4, 2, 1, 3 } }
};

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

    // Read and discard to newline
    static inline void readToNewline(ISstream& is)
    {
        char ch = '\n';
        do
        {
            is.get(ch);
        }
        while ((is) && ch != '\n');
    }

} // End namespace Foam


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

// Read two-line header
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_(BOUNDARY|CELL|VERTEX) [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  ...

\*---------------------------------------------------------------------------*/

bool Foam::fileFormats::STARCDCore::readHeader
(
    IFstream& is,
    const enum fileHeader header
)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    word magic;
    is >> magic;
    readToNewline(is);

    label majorVersion;
    is >> majorVersion;
    readToNewline(is);

    // Add other checks ...
    if (magic != fileHeaders_[header])
    {
        Info<< "Header mismatch " << fileHeaders_[header]
            << "  " << is.name()
            << nl;

        return false;
    }

    return true;
}


void Foam::fileFormats::STARCDCore::writeHeader
(
    Ostream& os,
    const enum fileHeader header
)
{
    os  << fileHeaders_[header] << nl
        << 4000
        << ' ' << 0
        << ' ' << 0
        << ' ' << 0
        << ' ' << 0
        << ' ' << 0
        << ' ' << 0
        << ' ' << 0
        << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileFormats::STARCDCore::starFileName
(
    const fileName& base,
    const enum fileExt ext
)
{
    return base + '.' + fileExtensions_[ext];
}


void Foam::fileFormats::STARCDCore::removeFiles(const fileName& base)
{
    Foam::rm(starFileName(base, VRT_FILE));
    Foam::rm(starFileName(base, CEL_FILE));
    Foam::rm(starFileName(base, BND_FILE));
    Foam::rm(starFileName(base, INP_FILE));
}


Foam::label Foam::fileFormats::STARCDCore::readPoints
(
    IFstream& is,
    List<point>& points,
    List<label>& ids
)
{
    label maxId = 0;
    token tok;

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    readHeader(is, HEADER_VRT);

    // Reuse memory if possible
    DynamicList<point> dynPoints(std::move(points));
    DynamicList<label> dynPointId(std::move(ids));  // STARCD index of points

    dynPoints.clear();
    dynPointId.clear();

    {
        scalar x, y, z;

        while (is.read(tok).good() && tok.isLabel())
        {
            const label starVertexId = tok.labelToken();

            is >> x >> y >> z;

            maxId = max(maxId, starVertexId);

            dynPoints.append(point(x, y, z));
            dynPointId.append(starVertexId);
        }
    }

    points.transfer(dynPoints);
    ids.transfer(dynPointId);

    return maxId;
}


void Foam::fileFormats::STARCDCore::writePoints
(
    Ostream& os,
    const UList<point>& points,
    const scalar scaleFactor
)
{
    writeHeader(os, HEADER_VRT);

    // Set the precision of the points data to 10
    os.precision(10);

    // Force decimal point for Fortran input
    os.setf(std::ios::showpoint);

    label starVertId = 1;  // 1-based vertex labels

    for (const point& p : points)
    {
        // Convert [m] -> [mm] etc
        os
            << starVertId << ' '
            << scaleFactor * p.x() << ' '
            << scaleFactor * p.y() << ' '
            << scaleFactor * p.z() << nl;

        ++starVertId;
    }

    os.flush();
}


// ************************************************************************* //
