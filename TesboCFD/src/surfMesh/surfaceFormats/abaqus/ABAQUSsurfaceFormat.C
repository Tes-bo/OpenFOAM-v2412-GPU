/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "ABAQUSsurfaceFormat.H"
#include "ListOps.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline Foam::label Foam::fileFormats::ABAQUSsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const Face& f,
    label elemId,
    const std::string& elsetName,
    bool header
)
{
    const label n = f.size();

    if (n == 4)
    {
        if (header)
        {
            os  << "*ELEMENT, TYPE=S4";

            if (!elsetName.empty())
            {
                os  << ", ELSET=" << elsetName;
            }
            os  << nl;
        }

        os  << "  "
            << (++elemId) << ','
            << (f[0] + 1) << ','
            << (f[1] + 1) << ','
            << (f[2] + 1) << ','
            << (f[3] + 1) << nl;
    }
    else
    {
        if (header)
        {
            os  << "*ELEMENT, TYPE=S3";

            if (!elsetName.empty())
            {
                os  << ", ELSET=" << elsetName;
            }
            os  << nl;
        }

        if (n == 3)
        {
            os  << "  "
                << (++elemId) << ','
                << (f[0] + 1) << ','
                << (f[1] + 1) << ','
                << (f[2] + 1) << nl;
        }
        else
        {
            // simple triangulation about f[0].
            // better triangulation should have been done before

            for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
            {
                const label fp2 = f.fcIndex(fp1);

                os  << "  "
                    << (++elemId) << ','
                    << (f[0] + 1) << ','
                    << (f[fp1] + 1) << ','
                    << (f[fp2] + 1) << nl;
            }
        }
    }

    return elemId;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::ABAQUSsurfaceFormat<Face>::ABAQUSsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::ABAQUSsurfaceFormat<Face>::read
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


    #ifdef FULLDEBUG
    ABAQUSCore::readHelper reader(true);  // Debugging verbosity
    #else
    ABAQUSCore::readHelper reader;
    #endif


    reader.read(is);

    // This needs more work

    // No solids
    reader.purge_solids();
    reader.compact_nodes();
    reader.renumber_elements_1to0();  // Convert 1-based -> 0-based


    // Convert connectivity to faces

    DynamicList<Face> dynFaces(reader.connectivity_.size());

    for (labelList& conn : reader.connectivity_)
    {
        dynFaces.append(Face(std::move(conn)));
    }


    // Rationalize the zones (elset)

    // Only retain element sets that are actually used
    labelHashSet elsetUsed(reader.elsetIds_);

    labelList newToOldZone(elsetUsed.sortedToc());

    // Extra safety
    if (newToOldZone.empty())
    {
        newToOldZone.push_back(0);
    }

    wordList  zoneNames(newToOldZone.size());
    labelList zoneSizes(newToOldZone.size(), Foam::zero{});

    Map<label> oldToNewZone(invertToMap(newToOldZone));

    forAllConstIters(reader.elsetMap_, iter)
    {
        const label zonei = oldToNewZone.lookup(iter.val(), -1);

        if (zonei >= 0)
        {
            zoneNames[zonei] = word::validate(iter.key());
        }
    }

    // No empty strings
    forAll(zoneNames, zonei)
    {
        if (zoneNames[zonei].empty())
        {
            zoneNames[zonei] = surfZoneIdentifier::defaultName(zonei);
        }
    }

    // Steal the elset Ids for our zones
    DynamicList<label> dynZones(std::move(reader.elsetIds_));

    // Renumber elset -> zoneId and increment the count
    for (label& zonei : dynZones)
    {
        zonei = oldToNewZone.lookup(zonei, 0);

        ++zoneSizes[zonei];
    }


    // Transfer to normal lists
    this->storedPoints().transfer(reader.points_);

    this->sortFacesAndStore
    (
        dynFaces,
        dynZones,
        reader.elemIds_,
        true // sorted
    );

    // Add zones (retaining empty ones)
    this->addZones(zoneSizes, zoneNames);
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::ABAQUSsurfaceFormat<Face>::write
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
    const UList<label>& elemIds  = surf.faceIds();

    // for no zones, suppress the group name
    const surfZoneList zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst, "")
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    // Possible to use faceIds?
    // - cannot if there are negative ids (eg, encoded solid/side)
    bool useOrigFaceIds =
    (
        !useFaceMap
     && elemIds.size() == faceLst.size()
     && !ListOps::found(elemIds, lessOp1<label>(0))
    );

    // Not possible with on-the-fly face decomposition
    if (useOrigFaceIds)
    {
        for (const auto& f : faceLst)
        {
            if (f.size() > 4)
            {
                useOrigFaceIds = false;
                break;
            }
        }
    }


    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }


    os  << "*HEADING" << nl;

    os  << nl
        << "**" << nl
        << "** Points" << nl
        << "**" << nl;

    writePoints(os, pointLst);

    os  << "**" << nl
        << "** Faces" << nl
        << "**" << nl
        << nl;


    // Simple tracking for change of element type/set
    labelPair prevOutput(-1, -1);

    label faceIndex = 0;
    label zoneIndex = 0;
    label elemId = 0;

    for (const surfZone& zone : zones)
    {
        for (label nLocal = zone.size(); nLocal--; ++faceIndex)
        {
            const label facei =
                (useFaceMap ? faceMap[faceIndex] : faceIndex);

            const Face& f = faceLst[facei];

            if (useOrigFaceIds)
            {
                elemId = elemIds[facei];
            }

            const label n = f.size();

            bool header =
                (prevOutput.first() != n || prevOutput.second() != zoneIndex);

            if (header)
            {
                // Update values
                prevOutput.first() = n;
                prevOutput.second() = zoneIndex;
            }

            elemId = writeShell(os, f, elemId, zone.name(), header);
        }

        ++zoneIndex;
    }

    os  << "**" << nl
        << "**" << nl;
}


// ************************************************************************* //
