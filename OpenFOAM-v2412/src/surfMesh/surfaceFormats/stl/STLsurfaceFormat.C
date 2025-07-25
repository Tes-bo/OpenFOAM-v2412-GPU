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

#include "STLsurfaceFormat.H"
#include "triangle.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::STLsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const UList<point>& pts,
    const Face& f
)
{
    // Calculate the normal ourselves, for flexibility and speed
    const vector norm =
        triPointRef::unitNormal(pts[f[0]], pts[f[1]], pts[f[2]]);

    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pts[f[0]];

    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        const label fp2 = f.fcIndex(fp1);

        // Write ASCII
        STLtriangle::write
        (
            os,
            norm,
            p0,
            pts[f[fp1]],
            pts[f[fp2]]
        );
    }
}


template<class Face>
inline void Foam::fileFormats::STLsurfaceFormat<Face>::writeShell
(
    ostream& os,
    const UList<point>& pts,
    const Face& f,
    const label zoneI
)
{
    // Calculate the normal ourselves, for flexibility and speed
    const vector norm =
        triPointRef::unitNormal(pts[f[0]], pts[f[1]], pts[f[2]]);

    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pts[f[0]];

    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        const label fp2 = f.fcIndex(fp1);

        // Write BINARY
        STLtriangle
        (
            norm,
            p0,
            pts[f[fp1]],
            pts[f[fp2]],
            zoneI
        ).write(os);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::STLsurfaceFormat<Face>::STLsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STLsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    // Clear everything
    this->clear();

    // Read in the values
    STLReader reader(filename);

    // Get the map for stitched surface points, with merge tolerance depending
    // on the input format
    labelList pointMap;
    const label nUniquePoints = reader.mergePointsMap(pointMap);

    const auto& readpts = reader.points();

    // Assign points
    pointField& pointLst = this->storedPoints();
    pointLst.setSize(nUniquePoints);
    forAll(readpts, pointi)
    {
        pointLst[pointMap[pointi]] = readpts[pointi];
    }

    // Retrieve the original zone information
    List<word>  names(std::move(reader.names()));
    List<label> sizes(std::move(reader.sizes()));
    List<label> zoneIds(std::move(reader.zoneIds()));

    // Generate the (sorted) faces
    List<Face> faceLst(zoneIds.size());

    if (reader.is_sorted())
    {
        // Already sorted - generate directly
        forAll(faceLst, facei)
        {
            const label startPt = 3*facei;
            faceLst[facei] = Face
            {
                pointMap[startPt],
                pointMap[startPt+1],
                pointMap[startPt+2]
            };
        }
    }
    else
    {
        // Determine the sorted order:
        // use sortedOrder directly (the intermediate list is discarded anyhow)
        labelList faceMap(sortedOrder(zoneIds));

        // Generate sorted faces
        forAll(faceMap, facei)
        {
            const label startPt = 3*faceMap[facei];
            faceLst[facei] = Face
            {
                pointMap[startPt],
                pointMap[startPt+1],
                pointMap[startPt+2]
            };
        }
    }
    zoneIds.clear();

    // Transfer
    this->storedFaces().transfer(faceLst);

    if (names.size())
    {
        this->addZones(sizes, names);
    }
    else
    {
        this->addZones(sizes);
    }
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeAscii
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    IOstreamOption::compressionType comp
)
{
    OFstream os(filename, IOstreamOption(IOstreamOption::ASCII, comp));
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    const UList<point>& pointLst = surf.points();
    const UList<Face>&   faceLst = surf.surfFaces();
    const UList<label>&  faceMap = surf.faceMap();

    const surfZoneList zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    label faceIndex = 0;
    for (const surfZone& zone : zones)
    {
        os << "solid " << zone.name() << nl;

        for (label nLocal = zone.size(); nLocal--; ++faceIndex)
        {
            const label facei =
                (useFaceMap ? faceMap[faceIndex] : faceIndex);

            const Face& f = faceLst[facei];

            writeShell(os, pointLst, f);
        }
        os << "endsolid " << zone.name() << endl;
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeBinary
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    std::ofstream os(filename, std::ios::binary);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    const UList<point>& pointLst = surf.points();
    const UList<Face>&   faceLst = surf.surfFaces();
    const UList<label>&  faceMap = surf.faceMap();

    const surfZoneList zones =
    (
        surf.surfZones().size() > 1
      ? surf.surfZones()
      : surfaceFormatsCore::oneZone(faceLst)
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    // Write the STL header
    unsigned int nTris = surf.nTriangles();
    STLCore::writeBinaryHeader(os, nTris);

    label faceIndex = 0;
    label zoneIndex = 0;
    for (const surfZone& zone : zones)
    {
        for (label nLocal = zone.size(); nLocal--; ++faceIndex)
        {
            const label facei =
                (useFaceMap ? faceMap[faceIndex] : faceIndex);

            const Face& f = faceLst[facei];

            writeShell(os, pointLst, f, zoneIndex);
        }

        ++zoneIndex;
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeAscii
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    IOstreamOption::compressionType comp
)
{
    const pointField& pointLst = surf.points();
    const UList<Face>& faceLst = surf.surfFaces();

    // A single zone - we can skip sorting
    if (surf.zoneToc().size() == 1)
    {
        OFstream os(filename, IOstreamOption(IOstreamOption::ASCII, comp));
        if (!os.good())
        {
            FatalErrorInFunction
                << "Cannot write file " << filename << nl
                << exit(FatalError);
        }

        os << "solid " << surf.zoneToc()[0].name() << nl;
        for (const Face& f : faceLst)
        {
            writeShell(os, pointLst, f);
        }
        os << "endsolid " << surf.zoneToc()[0].name() << nl;
    }
    else
    {
        labelList faceMap;
        List<surfZone> zoneLst = surf.sortedZones(faceMap);

        writeAscii
        (
            filename,
            MeshedSurfaceProxy<Face>
            (
                pointLst,
                faceLst,
                zoneLst,
                faceMap
            ),
            comp
        );
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeBinary
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    std::ofstream os(filename, std::ios::binary);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    const pointField& pointLst = surf.points();
    const UList<Face>& faceLst = surf.surfFaces();
    const UList<label>& zoneIds = surf.zoneIds();

    // Write the STL header
    unsigned int nTris = surf.nTriangles();
    STLCore::writeBinaryHeader(os, nTris);

    // Always write unsorted
    forAll(faceLst, facei)
    {
        writeShell
        (
            os,
            pointLst,
            faceLst[facei],
            zoneIds[facei]
        );
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary& options
)
{
    if
    (
        streamOpt.format() == IOstreamOption::BINARY
        // Detected "stlb" extension?
     || STLCore::isBinaryName(filename, STLCore::UNKNOWN)
    )
    {
        writeBinary(filename, surf);
    }
    else
    {
        writeAscii(filename, surf, streamOpt.compression());
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    const STLFormat format,
    IOstreamOption::compressionType comp
)
{
    if (STLCore::isBinaryName(filename, format))
    {
        writeBinary(filename, surf);
    }
    else
    {
        writeAscii(filename, surf, comp);
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary& options
)
{
    if
    (
        streamOpt.format() == IOstreamOption::BINARY
        // Detected "stlb" extension?
     || STLCore::isBinaryName(filename, STLCore::UNKNOWN)
    )
    {
        writeBinary(filename, surf);
    }
    else
    {
        writeAscii(filename, surf, streamOpt.compression());
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    const STLFormat format,
    IOstreamOption::compressionType comp
)
{
    if (STLCore::isBinaryName(filename, format))
    {
        writeBinary(filename, surf);
    }
    else
    {
        writeAscii(filename, surf, comp);
    }
}


// ************************************************************************* //
