/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "NASedgeFormat.H"
#include "IFstream.H"
#include "StringStream.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::NASedgeFormat::NASedgeFormat(const fileName& filename)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::NASedgeFormat::read
(
    const fileName& filename
)
{
    clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    DynamicList<point>  dynPoints;
    DynamicList<edge>   dynEdges;
    DynamicList<label>  pointId;     // Nastran index of points

    while (is.good())
    {
        string line;
        is.getLine(line);

        if (line.empty())
        {
            continue;  // Ignore empty
        }
        else if (line[0] == '$')
        {
            // Ignore comment
            continue;
        }

        // Check if character 72 is continuation
        if (line.size() > 72 && line[72] == '+')
        {
            line.resize(72);

            while (true)
            {
                string buf;
                is.getLine(buf);

                if (buf.size() > 72 && buf[72] == '+')
                {
                    line += buf.substr(8, 64);
                }
                else
                {
                    line += buf.substr(8);
                    break;
                }
            }
        }


        // Parsing position within current line
        std::string::size_type linei = 0;

        // Is free format if line contains a comma
        const bool freeFormat = line.contains(',');

        // First word (column 0-8)
        const word cmd(word::validate(nextNasField(line, linei, 8)));

        if (cmd == "CBEAM" || cmd == "CROD")
        {
            // Fixed format:
            //  8-16 : element id
            // 16-24 : group id
            // 24-32 : vertex
            // 32-40 : vertex

            // discard elementId
            (void) nextNasField(line, linei, 8, freeFormat);
            // discard groupId
            (void) nextNasField(line, linei, 8, freeFormat);

            label a = readLabel(nextNasField(line, linei, 8, freeFormat));
            label b = readLabel(nextNasField(line, linei, 8, freeFormat));

            dynEdges.emplace_back(a,b);
        }
        else if (cmd == "PLOTEL")
        {
            // Fixed format:
            //  8-16 : element id
            // 16-24 : vertex
            // 24-32 : vertex
            // 32-40 : vertex

            // discard elementId (8-16)
            (void) nextNasField(line, linei, 8, freeFormat);

            label a = readLabel(nextNasField(line, linei, 8, freeFormat));
            label b = readLabel(nextNasField(line, linei, 8, freeFormat));

            dynEdges.emplace_back(a,b);
        }
        else if (cmd == "GRID")
        {
            // Fixed (short) format:
            //  8-16 : point id
            // 16-24 : coordinate system (unsupported)
            // 24-32 : point x coordinate
            // 32-40 : point y coordinate
            // 40-48 : point z coordinate
            // 48-56 : displacement coordinate system (optional, unsupported)
            // 56-64 : single point constraints (optional, unsupported)
            // 64-70 : super-element id (optional, unsupported)

            label index = readLabel(nextNasField(line, linei, 8, freeFormat));
            (void) nextNasField(line, linei, 8, freeFormat);
            scalar x = readNasScalar(nextNasField(line, linei, 8, freeFormat));
            scalar y = readNasScalar(nextNasField(line, linei, 8, freeFormat));
            scalar z = readNasScalar(nextNasField(line, linei, 8, freeFormat));

            pointId.push_back(index);
            dynPoints.emplace_back(x, y, z);
        }
        else if (cmd == "GRID*")
        {
            // Long format is on two lines with '*' continuation symbol
            // on start of second line.
            // Typical line (spaces compacted)
            // GRID*      126   0 -5.55999875E+02 -5.68730474E+02
            // *         2.14897901E+02

            // Cannot be long format and free format at the same time!

            label index = readLabel(nextNasField(line, linei, 16)); // 8-24
            (void) nextNasField(line, linei, 16); // 24-40
            scalar x = readNasScalar(nextNasField(line, linei, 16)); // 40-56
            scalar y = readNasScalar(nextNasField(line, linei, 16)); // 56-72

            linei = 0; // restart at index 0
            is.getLine(line);
            if (line[0] != '*')
            {
                FatalErrorInFunction
                    << "Expected continuation symbol '*' when reading GRID*"
                    << " (double precision coordinate) format" << nl
                    << "Read:" << line << nl
                    << "File:" << is.name() << " line:" << is.lineNumber()
                    << exit(FatalError);
            }
            (void) nextNasField(line, linei, 8); // 0-8
            scalar z = readNasScalar(nextNasField(line, linei, 16)); // 8-16

            pointId.push_back(index);
            dynPoints.emplace_back(x, y, z);
        }
    }

    // transfer to normal lists
    storedPoints().transfer(dynPoints);

    // Build inverse mapping (NASTRAN pointId -> index)
    Map<label> mapPointId(invertToMap(pointId));

    // note which points were really used and which can be culled
    bitSet usedPoints(points().size());


    // Pass1: relabel edges
    // ~~~~~~~~~~~~~~~~~~~~
    for (edge& e : dynEdges)
    {
        e[0] = mapPointId[e[0]];
        e[1] = mapPointId[e[1]];

        usedPoints.set(e[0]);
        usedPoints.set(e[1]);
    }
    pointId.clearStorage();
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

            // map prev -> new id
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


// ************************************************************************* //
