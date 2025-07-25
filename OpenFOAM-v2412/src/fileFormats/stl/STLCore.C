/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "STLCore.H"
#include "OSspecific.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The number of bytes in the STL binary header
static constexpr const unsigned STLHeaderSize = 80;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Check if "SOLID" or "solid" appears as the first non-space content.
// Assume that any leading space is less than 75 chars or so, otherwise
// it is really bad input.
static bool startsWithSolid(const char header[STLHeaderSize])
{
    unsigned pos = 0;
    while (std::isspace(header[pos]) && pos < STLHeaderSize)
    {
        ++pos;
    }

    return
    (
        pos < (STLHeaderSize-5)  // At least 5 chars remaining
     && std::toupper(header[pos+0]) == 'S'
     && std::toupper(header[pos+1]) == 'O'
     && std::toupper(header[pos+2]) == 'L'
     && std::toupper(header[pos+3]) == 'I'
     && std::toupper(header[pos+4]) == 'D'
    );
}


// Check if file size appears to be reasonable for an STL binary file.
// Compare file size with that expected from number of tris
// If this is not sensible, it may be an ASCII file
//
// sizeof(STLtriangle) = 50 bytes [int16 + 4 * (3 float)]

inline static bool checkBinaryFileSize
(
    const int64_t nTris,
    const Foam::fileName& file
)
{
    // When checking the content size, account for the header size (80),
    // but ignore the nTris information (int32_t) to give some rounding

    const int64_t contentSize =
    (
        int64_t(Foam::fileSize(file))
      - int64_t(STLHeaderSize)
    );

    return
    (
        (contentSize >= 0)
     && (nTris >= contentSize/50)
     && (nTris <= contentSize/25)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::STLCore::isBinaryName
(
    const fileName& filename,
    const STLFormat format
)
{
    return
    (
        format == STLFormat::UNKNOWN
      ? filename.has_ext("stlb")
      : format == STLFormat::BINARY
    );
}


// Check binary by getting the header and number of facets
// this seems to work better than the old token-based method
// - using wordToken can cause an abort if non-word (binary) content
//   is detected ... this is not exactly what we want.
// - some programs (eg, PROSTAR) have 'solid' as the first word in
//   the binary header. This is just wrong and not our fault.
int Foam::fileFormats::STLCore::detectBinaryHeader
(
    const fileName& filename
)
{
    // Handle compressed (.gz) or uncompressed input files

    ifstreamPointer isPtr(filename);
    const bool unCompressed =
        (IOstreamOption::UNCOMPRESSED == isPtr.whichCompression());

    auto& is = *isPtr;

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << " or file " << filename + ".gz"
            << exit(FatalError);
    }

    // Read the STL header
    char header[STLHeaderSize];
    is.read(header, STLHeaderSize);

    // If the stream is bad, it can't be a binary STL
    if (!is.good() || startsWithSolid(header))
    {
        return 0;
    }


    // Read the number of triangles in the STL file
    // (note: read as signed int so we can check whether >2^31).
    //
    // With nTris == 2^31, file size is 107.37 GB !
    //
    // However, the limit is more likely caused by the number of points
    // that can be stored (label-size=32) when flattened for merging.
    // So more like 715.8M triangles (~35.8 GB)

    int32_t nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(int32_t));

    bool ok = (is && nTris >= 0);

    if (ok && unCompressed)
    {
        ok = checkBinaryFileSize(nTris, filename);
    }

    //if (ok)
    //{
    //    InfoErr<< "stlb : " << nTris << " triangles" << nl;
    //}

    // Return number of triangles if it appears to be BINARY and good.
    return (ok ? nTris : 0);
}


std::unique_ptr<std::istream>
Foam::fileFormats::STLCore::readBinaryHeader
(
    const fileName& filename,
    label& nTrisEstimated
)
{
    nTrisEstimated = 0;

    std::unique_ptr<std::istream> streamPtr;
    bool unCompressed(true);

    // Handle compressed (.gz) or uncompressed input files
    {
        ifstreamPointer isPtr(filename);
        unCompressed =
            (IOstreamOption::UNCOMPRESSED == isPtr.whichCompression());

        // Take ownership
        streamPtr.reset(isPtr.release());
    }
    auto& is = *streamPtr;

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << " or file " << filename + ".gz"
            << exit(FatalError);
    }


    // Read the STL header
    char header[STLHeaderSize];
    is.read(header, STLHeaderSize);

    // Check that stream is OK, if not this may be an ASCII file
    if (!is.good()) // could check again: startsWithSolid(header)
    {
        FatalErrorInFunction
            << "problem reading header, perhaps file is not binary "
            << exit(FatalError);
    }


    // Read the number of triangles in the STL file
    // (note: read as signed int so we can check whether >2^31).
    //
    // With nTris == 2^31, file size is 107.37 GB !
    //
    // However, the limit is more likely caused by the number of points
    // that can be stored (label-size=32) when flattened for merging.
    // So more like 715.8M triangles (~35.8 GB)

    int32_t nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(int32_t));

    bool ok = (is && nTris >= 0);

    if (ok && unCompressed)
    {
        ok = checkBinaryFileSize(nTris, filename);
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "problem reading number of triangles, perhaps file is not binary"
            << exit(FatalError);
    }

    nTrisEstimated = nTris;

    return streamPtr;
}


void Foam::fileFormats::STLCore::writeBinaryHeader
(
    ostream& os,
    uint32_t nTris
)
{
    // STL header with extra information about nTris
    char header[STLHeaderSize];
    ::snprintf(header, STLHeaderSize, "STL binary file %u facets", nTris);

    // Fill trailing with zeroes (to avoid writing junk)
    for (size_t i = strlen(header); i < STLHeaderSize; ++i)
    {
        header[i] = 0;
    }

    os.write(header, STLHeaderSize);
    os.write(reinterpret_cast<char*>(&nTris), sizeof(uint32_t));
}


// ************************************************************************* //
