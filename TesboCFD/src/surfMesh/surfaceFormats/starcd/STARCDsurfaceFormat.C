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

#include "STARCDsurfaceFormat.H"
#include "ListOps.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::STARCDsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const Face& f,
    const label cellId,
    const label cellTableId
)
{
    os  << (cellId + 1)
        << ' ' << starcdShell       // 3(shell) shape
        << ' ' << f.size()
        << ' ' << (cellTableId + 1)
        << ' ' << starcdShellType;  // 4(shell)

    // Primitives have <= 8 vertices, but prevent overrun anyhow
    // indent following lines for ease of reading
    label count = 0;
    for (const label pointi : f)
    {
        if ((count % 8) == 0)
        {
            os  << nl << "  " << (cellId + 1);
        }
        os  << ' ' << (pointi + 1);
        ++count;
    }
    os  << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::STARCDsurfaceFormat<Face>::STARCDsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STARCDsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    // Clear everything
    this->clear();

    const fileName prefix(filename.lessExt());

    // Read cellTable names (if possible)
    Map<word> cellTableLookup = readInpCellTable
    (
        IFstream(starFileName(prefix, STARCDCore::INP_FILE))()
    );


    // STARCD index of points
    List<label> pointId;

    // read points from .vrt file
    readPoints
    (
        IFstream(starFileName(prefix, STARCDCore::VRT_FILE))(),
        this->storedPoints(),
        pointId
    );

    // Build inverse mapping (STARCD pointId -> index)
    Map<label> mapPointId(invertToMap(pointId));
    pointId.clear();


    // Read .cel file
    // ~~~~~~~~~~~~~~
    IFstream is(starFileName(prefix, STARCDCore::CEL_FILE));
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name() << nl
            << exit(FatalError);
    }

    readHeader(is, STARCDCore::HEADER_CEL);

    DynamicList<label> dynElemId;  // STARCD element id (1-based)
    DynamicList<Face>  dynFaces;

    DynamicList<label> dynZones;
    DynamicList<word>  dynNames;
    DynamicList<label> dynSizes;
    Map<label> lookup;

    // Assume the cellTableIds are not intermixed
    bool sorted = true;
    label zoneId = 0;

    // Element id gets trashed with decompose into a triangle!
    bool ignoreElemId = false;

    label ignoredLabel, shapeId, nLabels, cellTableId, typeId;
    DynamicList<label> vertexLabels(64);

    token tok;

    while (is.read(tok).good() && tok.isLabel())
    {
        // First token is the element id (1-based)
        label elemId = tok.labelToken();

        is  >> shapeId
            >> nLabels
            >> cellTableId
            >> typeId;

        vertexLabels.clear();
        vertexLabels.reserve(nLabels);

        // Read indices - max 8 per line
        for (label i = 0; i < nLabels; ++i)
        {
            label vrtId;
            if ((i % 8) == 0)
            {
                is >> ignoredLabel; // Skip cellId for continuation lines
            }
            is >> vrtId;

            // Convert original vertex id to point label
            vertexLabels.append(mapPointId[vrtId]);
        }

        if (typeId == starcdShellType)
        {
            // Convert cellTableId to zoneId
            const auto iterGroup = lookup.cfind(cellTableId);
            if (iterGroup.good())
            {
                if (zoneId != *iterGroup)
                {
                    // cellTableIds are intermixed
                    sorted = false;
                }
                zoneId = *iterGroup;
            }
            else
            {
                zoneId = dynSizes.size();
                lookup.insert(cellTableId, zoneId);

                const auto iterTableName = cellTableLookup.cfind(cellTableId);

                if (iterTableName.good())
                {
                    dynNames.append(*iterTableName);
                }
                else
                {
                    dynNames.append("cellTable_" + ::Foam::name(cellTableId));
                }

                dynSizes.append(0);
            }

            SubList<label> vertices(vertexLabels, vertexLabels.size());
            if (faceTraits<Face>::isTri() && nLabels > 3)
            {
                // The face needs triangulation
                ignoreElemId = true;
                dynElemId.clear();

                face f(vertices);

                faceList trias(f.nTriangles());
                label nTri = 0;
                f.triangles(this->points(), nTri, trias);

                for (const face& tri : trias)
                {
                    // A triangular 'face', convert to 'triFace' etc
                    dynFaces.append(Face(tri));
                    dynZones.append(zoneId);
                    dynSizes[zoneId]++;
                }
            }
            else if (nLabels >= 3)
            {
                --elemId;   // Convert 1-based -> 0-based
                dynElemId.append(elemId);

                dynFaces.append(Face(vertices));
                dynZones.append(zoneId);
                dynSizes[zoneId]++;
            }
        }
    }
    mapPointId.clear();


    if (ignoreElemId)
    {
        dynElemId.clear();
    }


    this->sortFacesAndStore(dynFaces, dynZones, dynElemId, sorted);

    // Add zones (retaining empty ones)
    this->addZones(dynSizes, dynNames);
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::STARCDsurfaceFormat<Face>::write
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

    const surfZoneList zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    // Possible to use faceIds?
    // - cannot if there are negative ids (eg, encoded solid/side)
    const bool useOrigFaceIds =
    (
        !useFaceMap
     && elemIds.size() == faceLst.size()
     && !ListOps::found(elemIds, lessOp1<label>(0))
    );


    const fileName prefix(filename.lessExt());

    // The .vrt file
    {
        OFstream os(starFileName(prefix, STARCDCore::VRT_FILE), streamOpt);
        writePoints(os, pointLst);
    }

    // The .cel file
    OFstream os(starFileName(prefix, STARCDCore::CEL_FILE), streamOpt);
    writeHeader(os, STARCDCore::HEADER_CEL);

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

            writeShell(os, f, elemId, zoneIndex);
            ++elemId;
        }

        ++zoneIndex;
    }

    // Simple .inp file - always UNCOMPRESSED
    {
        OFstream os(starFileName(prefix, STARCDCore::INP_FILE));

        writeCase
        (
            os,
            pointLst,
            faceLst.size(),
            zones
        );
    }
}


// ************************************************************************* //
