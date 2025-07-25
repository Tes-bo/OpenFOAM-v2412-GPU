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

#include "TRIsurfaceFormat.H"
#include "TRIReader.H"
#include "OFstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::TRIsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const UList<point>& pts,
    const Face& f,
    const label zoneI
)
{
    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pts[f[0]];
    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        const label fp2 = f.fcIndex(fp1);

        const point& p1 = pts[f[fp1]];
        const point& p2 = pts[f[fp2]];

        os  << p0.x() << ' ' << p0.y() << ' ' << p0.z() << ' '
            << p1.x() << ' ' << p1.y() << ' ' << p1.z() << ' '
            << p2.x() << ' ' << p2.y() << ' ' << p2.z() << ' '
            // zone as colour
            << "0x" << hex << zoneI << dec << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::TRIsurfaceFormat<Face>::TRIsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::TRIsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    // Clear everything
    this->clear();

    // Read in the values
    TRIReader reader(filename);

    // Get the map for stitched surface points
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
        // use sortedOrder directly (the intermediate list is discared anyhow)
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

    this->addZones(sizes);
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::TRIsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary&
)
{
    // ASCII only, allow output compression
    streamOpt.format(IOstreamOption::ASCII);

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

    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

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
void Foam::fileFormats::TRIsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary&
)
{
    // ASCII only, allow output compression
    streamOpt.format(IOstreamOption::ASCII);

    const UList<point>& pointLst = surf.points();
    const UList<Face>&   faceLst = surf.surfFaces();

    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    // A single zone needs no sorting
    if (surf.zoneToc().size() == 1)
    {
        const UList<label>& zoneIds = surf.zoneIds();

        forAll(faceLst, facei)
        {
            writeShell(os, pointLst, faceLst[facei], zoneIds[facei]);
        }
    }
    else
    {
        labelList faceMap;
        List<surfZone> zoneLst = surf.sortedZones(faceMap);

        label faceIndex = 0;
        label zoneIndex = 0;

        for (const surfZone& zone : zoneLst)
        {
            for (label nLocal = zone.size(); nLocal--; ++faceIndex)
            {
                const label facei = faceMap[faceIndex];

                const Face& f = faceLst[facei];

                writeShell(os, pointLst, f, zoneIndex);
            }

            ++zoneIndex;
        }
    }
}


// ************************************************************************* //
