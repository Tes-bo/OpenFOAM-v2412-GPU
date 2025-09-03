/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "STARCDedgeFormat.H"
#include "ListOps.H"
#include "clock.H"
#include "bitSet.H"
#include "StringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::fileFormats::STARCDedgeFormat::writeLines
(
    Ostream& os,
    const edgeList& edges,
    label starCellId
)
{
    starCellId = max(1, starCellId);   // Enforce 1-based cellId

    for (const edge& e : edges)
    {
        os  << starCellId
            << ' ' << starcdLine         // 2(line) shape
            << ' ' << e.size()
            << ' ' << 401                // arbitrary value
            << ' ' << starcdLineType;    // 5(line)

        os  << nl
            << "  " << starCellId << "  "
            << (e[0]+1) << "  " << (e[1]+1) << nl;

        ++starCellId;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::STARCDedgeFormat::writeCase
(
    Ostream& os,
    const pointField& pointLst,
    const label nEdges
)
{
    const word caseName = os.name().stem();

    os  << "! STARCD file written " << clock::dateTime().c_str() << nl
        << "! " << pointLst.size() << " points, " << nEdges << " lines" << nl
        << "! case " << caseName << nl
        << "! ------------------------------" << nl;

//     forAll(zoneLst, zoneI)
//     {
//         os  << "ctable " << zoneI + 1 << " line" << nl
//             << "ctname " << zoneI + 1 << " "
//             << zoneLst[zoneI].name() << nl;
//     }

    os  << "! ------------------------------" << nl
        << "*set icvo mxv - 1" << nl
        << "vread " << caseName << ".vrt icvo,,,coded" << nl
        << "cread " << caseName << ".cel icvo,,,add,coded" << nl
        << "*set icvo" << nl
        << "! end" << nl;

    os.flush();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::STARCDedgeFormat::STARCDedgeFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::STARCDedgeFormat::read
(
    const fileName& filename
)
{
    clear();

    const fileName prefix(filename.lessExt());

    // STARCD index of points
    List<label> pointId;

    // Read points from .vrt file
    readPoints
    (
        IFstream(starFileName(prefix, STARCDCore::VRT_FILE))(),
        storedPoints(),
        pointId
    );

    // Build inverse mapping (STARCD pointId -> index)
    Map<label> mapPointId(invertToMap(pointId));
    pointId.clear();

    // Note which points were really used and which can be culled
    bitSet usedPoints(points().size());


    // Read .cel file
    // ~~~~~~~~~~~~~~
    IFstream is(starFileName(prefix, STARCDCore::CEL_FILE));
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    readHeader(is, STARCDCore::HEADER_CEL);

    DynamicList<edge> dynEdges;

    label ignoredLabel, shapeId, nLabels, cellTableId, typeId;
    DynamicList<label> vertexLabels(64);

    token tok;

    while (is.read(tok).good() && tok.isLabel())
    {
        // const label starCellId = tok.labelToken();
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

        if (typeId == starcdLineType)
        {
            if (vertexLabels.size() >= 2)
            {
                dynEdges.append(edge(vertexLabels[0], vertexLabels[1]));

                usedPoints.set(vertexLabels[0]);
                usedPoints.set(vertexLabels[1]);
            }
        }
    }

    mapPointId.clear();

    // Not all points were used, subset/cull them accordingly
    if (!usedPoints.all())
    {
        label nUsed = 0;

        pointField& pts = storedPoints();
        for (const label pointi : usedPoints)
        {
            if (nUsed != pointi)
            {
                pts[nUsed] = pts[pointi];
            }

            // Map prev -> new id
            mapPointId.set(pointi, nUsed);

            ++nUsed;
        }
        pts.resize(nUsed);

        // Renumber edge vertices
        for (edge& e : dynEdges)
        {
            e[0] = mapPointId[e[0]];
            e[1] = mapPointId[e[1]];
        }
    }

    storedEdges().transfer(dynEdges);

    return true;
}


void Foam::fileFormats::STARCDedgeFormat::write
(
    const fileName& filename,
    const edgeMesh& mesh,
    IOstreamOption streamOpt,
    const dictionary&
)
{
    // ASCII only, allow output compression
    streamOpt.format(IOstreamOption::ASCII);

    const pointField& pointLst = mesh.points();
    const edgeList& edgeLst = mesh.edges();

    const fileName prefix(filename.lessExt());

    // The .vrt file
    {
        OFstream os(starFileName(prefix, STARCDCore::VRT_FILE), streamOpt);
        writePoints(os, pointLst);
    }

    // The .cel file
    {
        OFstream os(starFileName(prefix, STARCDCore::CEL_FILE), streamOpt);
        writeHeader(os, STARCDCore::HEADER_CEL);
        writeLines(os, edgeLst);
    }

    // Write a simple .inp file. Never compressed
    writeCase
    (
        OFstream(starFileName(prefix, STARCDCore::INP_FILE))(),
        pointLst,
        edgeLst.size()
    );
}


// ************************************************************************* //
