/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "STLReader.H"
#include "STLAsciiParse.H"
#include "Map.H"
#include "IFstream.H"
#include "mergePoints.H"

#undef DEBUG_STLBINARY

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

int Foam::fileFormats::STLReader::parserType
(
    Foam::debug::optimisationSwitch("fileFormats::stl", 0)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::STLReader::transfer
(
    Detail::STLAsciiParse& parsed
)
{
    sorted_ = parsed.is_sorted();

    points_.transfer(parsed.points());
    zoneIds_.transfer(parsed.facets());
    names_.transfer(parsed.names());
    sizes_.transfer(parsed.sizes());

    format_ = STLFormat::ASCII;

    parsed.clear();
}


bool Foam::fileFormats::STLReader::readASCII
(
    const fileName& filename
)
{
    // No runtime selection of parser (only via optimisationSwitch)
    // this is something that is infrequently changed.
    if (parserType == 1)
    {
        return readAsciiRagel(filename);
    }
    else if (parserType == 2)
    {
        return readAsciiManual(filename);
    }
    return readAsciiFlex(filename);
}


bool Foam::fileFormats::STLReader::readBINARY
(
    const fileName& filename
)
{
    sorted_ = true;
    format_ = STLFormat::UNKNOWN;

    label nTris = 0;
    std::unique_ptr<std::istream> streamPtr
    {
        readBinaryHeader(filename, nTris)
    };

    if (!streamPtr)
    {
        FatalErrorInFunction
            << "Error reading file " << filename
            << " or file " << filename + ".gz"
            << exit(FatalError);
    }
    auto& is = *streamPtr;

#ifdef DEBUG_STLBINARY
    Info<< "# " << nTris << " facets" << endl;
    label prevZone = -1;
#endif

    points_.setSize(3*nTris);
    zoneIds_.setSize(nTris);

    Map<label> lookup;
    DynamicList<label> dynSizes;

    label ptI = 0;
    label zoneI = -1;
    forAll(zoneIds_, facei)
    {
        // Read STL triangle
        STLtriangle stlTri(is);

        // transcribe the vertices of the STL triangle -> points
        points_[ptI++] = stlTri.a();
        points_[ptI++] = stlTri.b();
        points_[ptI++] = stlTri.c();

        // interpret STL attribute as a zone
        const label origId = stlTri.attrib();

        auto fnd = lookup.cfind(origId);
        if (fnd.good())
        {
            if (zoneI != *fnd)
            {
                sorted_ = false; // Group appeared out of order
            }
            zoneI = *fnd;
        }
        else
        {
            zoneI = dynSizes.size();
            lookup.insert(origId, zoneI);
            dynSizes.push_back(0);
        }

        zoneIds_[facei] = zoneI;
        dynSizes[zoneI]++;

#ifdef DEBUG_STLBINARY
        if (prevZone != zoneI)
        {
            if (prevZone != -1)
            {
                Info<< "endsolid zone" << prevZone << nl;
            }
            prevZone = zoneI;

            Info<< "solid zone" << prevZone << nl;
        }

        stlTri.print(Info);
#endif
    }

#ifdef DEBUG_STLBINARY
    if (prevZone != -1)
    {
        Info<< "endsolid zone" << prevZone << nl;
    }
#endif

    names_.clear();
    sizes_.transfer(dynSizes);

    format_ = STLFormat::BINARY;
    return true;
}


bool Foam::fileFormats::STLReader::readFile
(
    const fileName& filename,
    const STLFormat format
)
{
    if
    (
        format == STLFormat::UNKNOWN
      ? detectBinaryHeader(filename)
      : format == STLFormat::BINARY
    )
    {
        return readBINARY(filename);
    }
    else
    {
        return readASCII(filename);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::STLReader::STLReader
(
    const fileName& filename
)
:
    sorted_(true),
    points_(),
    zoneIds_(),
    names_(),
    sizes_(),
    format_(STLFormat::UNKNOWN)
{
    // Auto-detect ASCII/BINARY format
    readFile(filename, STLFormat::UNKNOWN);
}


Foam::fileFormats::STLReader::STLReader
(
    const fileName& filename,
    const STLFormat format
)
:
    sorted_(true),
    points_(),
    zoneIds_(),
    names_(),
    sizes_(),
    format_(STLFormat::UNKNOWN)
{
    // Manually specified ASCII/BINARY format
    readFile(filename, format);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fileFormats::STLReader::clear()
{
    sorted_ = true;
    points_.clear();
    zoneIds_.clear();
    names_.clear();
    sizes_.clear();
    format_ = STLFormat::UNKNOWN;
}


Foam::label Foam::fileFormats::STLReader::mergePointsMap
(
    labelList& pointMap
) const
{
    // With the merge distance depending on the input format (ASCII | BINARY),
    // but must be independent of WM_SP or WM_DP flag.
    // - floatScalarSMALL  = 1e-6
    // - doubleScalarSMALL = 1e-15

    return mergePointsMap
    (
        (format_ == STLFormat::BINARY ? 10 : 100) * doubleScalarSMALL,
        pointMap
    );
}


Foam::label Foam::fileFormats::STLReader::mergePointsMap
(
    const scalar mergeTol,
    labelList& pointMap
) const
{
    return Foam::mergePoints
    (
        points_,
        mergeTol,
        false, // verbose
        pointMap
    );
}


// ************************************************************************* //
