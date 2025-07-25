/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "OBJsurfaceFormat.H"
#include "clock.H"
#include "Fstream.H"
#include "stringOps.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::OBJsurfaceFormat<Face>::OBJsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::OBJsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    // Clear everything
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename << nl
            << exit(FatalError);
    }

    // Assume that the groups are not intermixed
    // Place faces without a group in zone0.
    // zoneId = -1 to signal uninitialized
    label zoneId = -1;
    bool sorted = true;

    DynamicList<point> dynPoints;
    DynamicList<label> dynVerts;

    DynamicList<label> dynElemId; // unused
    DynamicList<Face>  dynFaces;

    DynamicList<word>  dynNames;
    DynamicList<label> dynZones;
    DynamicList<label> dynSizes;

    HashTable<label>   groupLookup;

    while (is.good())
    {
        string line = this->getLineNoComment(is);

        // Line continuations
        while (line.removeEnd('\\'))
        {
            line += this->getLineNoComment(is);
        }

        const auto tokens = stringOps::splitSpace(line);

        // Require command and some arguments
        if (tokens.size() < 2)
        {
            continue;
        }

        const word cmd = word::validate(tokens[0]);

        if (cmd == "v")
        {
            // Vertex
            // v x y z

            dynPoints.append
            (
                point
                (
                    readScalar(tokens[1]),
                    readScalar(tokens[2]),
                    readScalar(tokens[3])
                )
            );
        }
        else if (cmd == "g")
        {
            // Grouping
            // g name

            const word groupName = word::validate(tokens[1]);
            const auto iterGroup = groupLookup.cfind(groupName);

            if (iterGroup.good())
            {
                if (zoneId != *iterGroup)
                {
                    sorted = false; // Group appeared out of order
                }
                zoneId = *iterGroup;
            }
            else
            {
                zoneId = dynSizes.size();
                groupLookup.insert(groupName, zoneId);
                dynNames.append(groupName);
                dynSizes.append(0);
            }
        }
        else if (cmd == "f")
        {
            // Face
            // f v1 v2 v3 ...
            // OR
            // f v1/vt1 v2/vt2 v3/vt3 ...

            // Ensure it has as valid grouping
            if (zoneId < 0)
            {
                zoneId = 0;
                groupLookup.insert("zone0", 0);
                dynNames.append("zone0");
                dynSizes.append(0);
            }

            dynVerts.clear();

            bool first = true;
            for (const auto& tok : tokens)
            {
                if (first)
                {
                    // skip initial "f" or "l"
                    first = false;
                    continue;
                }

                const string vrtSpec(tok);
                const auto slash = vrtSpec.find('/');

                const label vertId =
                (
                    slash != string::npos
                  ? readLabel(vrtSpec.substr(0, slash))
                  : readLabel(vrtSpec)
                );

                dynVerts.append(vertId - 1);
            }

            const labelUList& f = dynVerts;

            if (faceTraits<Face>::isTri() && f.size() > 3)
            {
                // simple face triangulation about f[0]
                // points may be incomplete
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    const label fp2 = f.fcIndex(fp1);

                    dynFaces.append(Face{f[0], f[fp1], f[fp2]});
                    dynZones.append(zoneId);
                    dynSizes[zoneId]++;
                }
            }
            else if (f.size() >= 3)
            {
                dynFaces.append(Face(f));
                dynZones.append(zoneId);
                dynSizes[zoneId]++;
            }
        }
    }


    // Transfer to normal lists
    this->storedPoints().transfer(dynPoints);

    this->sortFacesAndStore(dynFaces, dynZones, dynElemId, sorted);

    // Add zones (retaining empty ones)
    this->addZones(dynSizes, dynNames);
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::OBJsurfaceFormat<Face>::write
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
    const UList<Face>&  faceLst  = surf.surfFaces();
    const UList<label>& faceMap  = surf.faceMap();

    // for no zones, suppress the group name
    const surfZoneList zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst, "")
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


    os  << "# Wavefront OBJ file written " << clock::dateTime().c_str() << nl
        << "o " << os.name().stem() << nl
        << nl
        << "# points : " << pointLst.size() << nl
        << "# faces  : " << faceLst.size() << nl
        << "# zones  : " << zones.size() << nl;

    // Print zone names as comment
    forAll(zones, zonei)
    {
        os  << "#   " << zonei << "  " << zones[zonei].name()
            << "  (nFaces: " << zones[zonei].size() << ")" << nl;
    }

    os  << nl
        << "# <points count=\"" << pointLst.size() << "\">" << nl;

    // Write vertex coords
    for (const point& pt : pointLst)
    {
        os  << "v " << pt.x() << ' '  << pt.y() << ' '  << pt.z() << nl;
    }

    os  << "# </points>" << nl
        << nl
        << "# <faces count=\"" << faceLst.size() << "\">" << nl;


    label faceIndex = 0;

    for (const surfZone& zone : zones)
    {
        if (zone.name().size())
        {
            os << "g " << zone.name() << nl;
        }

        for (label nLocal = zone.size(); nLocal--; ++faceIndex)
        {
            const label facei =
                (useFaceMap ? faceMap[faceIndex] : faceIndex);

            const Face& f = faceLst[facei];

            os << 'f';
            for (const label verti : f)
            {
                os << ' ' << (verti + 1);
            }
            os << nl;
        }
    }
    os << "# </faces>" << nl;
}


// ************************************************************************* //
