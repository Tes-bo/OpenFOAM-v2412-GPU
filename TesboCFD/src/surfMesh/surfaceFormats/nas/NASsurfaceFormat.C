/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Notes (Reader)

    Nastran does not directly support any names, but ANSA and Hypermesh
    have different ways to get around that problem by using/misusing
    comment lines.

Hypermesh extension (last verified approx. 2010)

    $HMNAME COMP                   1"some-part-name"

ANSA extension (legacy)

    line 1: $ANSA_NAME;<int>;<word>;
    line 2: $some-part-name

ANSA extension (19.0.1)

    line 1: $ANSA_NAME;<int>;PSHELL;~
    line 2: $some-part-name

    These seem to appear immediately before the corrsponding PSHELL

ANSA extension (23.1.0)

    $ANSA_NAME_COMMENT;<int>;PSHELL;some-part-name;; ...something trailing...

    These seem to appear as footer data, but could presumably appear anywhere.

Random extension (not sure where this arises)

    $some-part-name
    PSHELL  203101  1

    These seemingly random comments appear immediately before the PSHELL entry.

\*---------------------------------------------------------------------------*/

#include "NASsurfaceFormat.H"
#include "ListOps.H"
#include "Fstream.H"
#include "IOmanip.H"
#include "faceTraits.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline Foam::label Foam::fileFormats::NASsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const Face& f,
    label elemId,
    const label groupId
)
{
    const label n = f.size();

    if (n == 3)
    {
        os  << "CTRIA3" << ','
            << (++elemId) << ','
            << (groupId + 1) << ','
            << (f[0] + 1) << ','
            << (f[1] + 1) << ','
            << (f[2] + 1) << nl;
    }
    else if (n == 4)
    {
        os  << "CQUAD4" << ','
            << (++elemId) << ','
            << (groupId + 1) << ','
            << (f[0] + 1) << ','
            << (f[1] + 1) << ','
            << (f[2] + 1) << ','
            << (f[3] + 1) << nl;
    }
    else
    {
        // simple triangulation about f[0].
        // better triangulation should have been done before
        for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
        {
            const label fp2 = f.fcIndex(fp1);

            os  << "CTRIA3" << ','
                << (++elemId) << ','
                << (groupId + 1) << ','
                << (f[0] + 1) << ','
                << (f[fp1] + 1) << ','
                << (f[fp2] + 1) << nl;
        }
    }

    return elemId;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::NASsurfaceFormat<Face>::NASsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::NASsurfaceFormat<Face>::read
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

    DynamicList<label>  pointId;    // Nastran point id (1-based)
    DynamicList<point>  dynPoints;

    DynamicList<label>  dynElemId;  // Nastran element id (1-based)
    DynamicList<Face>   dynFaces;
    DynamicList<label>  dynZones;
    DynamicList<label>  dynSizes;

    Map<label>          zoneLookup;

    // Assume that the groups are not intermixed
    label zoneId = 0;
    bool sorted = true;

    // Element id gets trashed with decompose into a triangle!
    bool ignoreElemId = false;

    // Name for face group (limited to PSHELL)
    Map<word> nameLookup;

    // A single warning per unrecognized command
    wordHashSet unhandledCmd;

    // The line to parse
    string line;

    // The last comment line seen (immediately before a 'real' command)
    string lastComment;

    while (is.good())
    {
        is.getLine(line);

        if (NASCore::debug > 1) Info<< "Process: " << line << nl;

        // ANSA extension(s)
        if (line.starts_with("$ANSA_NAME"))
        {
            // Keep empty elements when splitting
            const auto args =
                stringOps::split<std::string>(line, ';', 0, true);

            if (args.size() > 4 && line.starts_with("$ANSA_NAME_COMMENT"))
            {
                // This type of content
                // $ANSA_NAME_COMMENT;93000;PSHELL;SLIP;;NO;NO;NO;NO;

                label groupId = 0;

                if (readLabel(args[1], groupId) && (args[2] == "PSHELL"))
                {
                    word groupName = word::validate(args[3]);

                    if (!groupName.empty())
                    {
                        DebugInfo
                            << "PSHELL:" << groupId
                            << " = " << groupName << nl;

                        nameLookup.emplace(groupId, std::move(groupName));
                    }
                }

                // Handled (or ignored)
                continue;
            }
            else if (args.size() >= 3 && (args[0] == "$ANSA_NAME"))
            {
                // This type of content

                // line 1: $ANSA_NAME;<int>;PSHELL;~
                // line 2: $some-part-name

                label groupId = 0;

                if (readLabel(args[1], groupId) && (args[2] == "PSHELL"))
                {
                    // Fetch the next line
                    is.getLine(line);
                    line.removeEnd('\r');  // Possible CR-NL

                    word groupName;
                    if (line.starts_with('$'))
                    {
                        groupName = word::validate(line.substr(1));
                    }

                    if (!groupName.empty())
                    {
                        DebugInfo
                            << "PSHELL:" << groupId
                            << " = " << groupName << nl;

                        nameLookup.emplace(groupId, std::move(groupName));
                    }
                }
            }

            // Drop through in case the second line read was not a comment !
        }
        else if (line.starts_with("$HMNAME COMP"))
        {
            // HYPERMESH extension
            // This type of content
            // $HMNAME COMP                   1"partName"
            // [NB: first entry is fixed record length of 32]

            auto dquote = line.find('"', 12);  // Beyond '$HMNAME COMP'

            label groupId = 0;

            if
            (
                dquote != std::string::npos
             && readLabel(line.substr(12, (dquote - 12)), groupId)
            )
            {
                // word::validate automatically removes quotes too
                word groupName = word::validate(line.substr(dquote));

                if (!groupName.empty())
                {
                    DebugInfo
                        << "HMNAME group " << groupId
                        << " => " << groupName << nl;

                    nameLookup.emplace(groupId, std::move(groupName));
                }
            }

            continue;  // Handled
        }

        if (line.empty())
        {
            continue;  // Ignore empty
        }
        else if (line[0] == '$')
        {
            // Retain comment (see notes above about weird formats...)
            lastComment = line;
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

        if (cmd == "CTRIA3")
        {
            // Fixed format:
            //  8-16 : element id
            // 16-24 : group id
            // 24-32 : vertex
            // 32-40 : vertex
            // 40-48 : vertex

            label elemId = readLabel(nextNasField(line, linei, 8, freeFormat));
            label groupId = readLabel(nextNasField(line, linei, 8, freeFormat));
            const auto a = readLabel(nextNasField(line, linei, 8, freeFormat));
            const auto b = readLabel(nextNasField(line, linei, 8, freeFormat));
            const auto c = readLabel(nextNasField(line, linei, 8, freeFormat));

            // Convert groupId into zoneId
            const auto iterZone = zoneLookup.cfind(groupId);
            if (iterZone.good())
            {
                if (zoneId != iterZone.val())
                {
                    // PSHELL types are intermixed
                    sorted = false;
                }
                zoneId = iterZone.val();
            }
            else
            {
                zoneId = dynSizes.size();
                zoneLookup.insert(groupId, zoneId);
                dynSizes.push_back(0);
                // Info<< "zone" << zoneId << " => group " << groupId <<nl;
            }

            --elemId;   // Convert 1-based -> 0-based
            dynElemId.push_back(elemId);
            dynFaces.push_back(Face{a, b, c});
            dynZones.push_back(zoneId);
            dynSizes[zoneId]++;
        }
        else if (cmd == "CQUAD4")
        {
            // Fixed format:
            //  8-16 : element id
            // 16-24 : group id
            // 24-32 : vertex
            // 32-40 : vertex
            // 40-48 : vertex
            // 48-56 : vertex

            label elemId = readLabel(nextNasField(line, linei, 8, freeFormat));
            label groupId = readLabel(nextNasField(line, linei, 8, freeFormat));
            const auto a = readLabel(nextNasField(line, linei, 8, freeFormat));
            const auto b = readLabel(nextNasField(line, linei, 8, freeFormat));
            const auto c = readLabel(nextNasField(line, linei, 8, freeFormat));
            const auto d = readLabel(nextNasField(line, linei, 8, freeFormat));

            // Convert groupId into zoneId
            const auto iterZone = zoneLookup.cfind(groupId);
            if (iterZone.good())
            {
                if (zoneId != iterZone.val())
                {
                    // PSHELL types are intermixed
                    sorted = false;
                }
                zoneId = iterZone.val();
            }
            else
            {
                zoneId = dynSizes.size();
                zoneLookup.insert(groupId, zoneId);
                dynSizes.push_back(0);
                // Info<< "zone" << zoneId << " => group " << groupId <<nl;
            }

            if (faceTraits<Face>::isTri())
            {
                ignoreElemId = true;
                dynElemId.clear();

                dynFaces.push_back(Face{a, b, c});
                dynFaces.push_back(Face{c, d, a});
                dynZones.push_back(zoneId);
                dynZones.push_back(zoneId);
                dynSizes[zoneId] += 2;
            }
            else
            {
                --elemId;   // Convert 1-based -> 0-based

                dynElemId.push_back(elemId);
                dynFaces.push_back(Face{a,b,c,d});
                dynZones.push_back(zoneId);
                dynSizes[zoneId]++;
            }
        }
        else if (cmd == "GRID")
        {
            // Fixed (short) format:
            //  8-16 : point id
            // 16-24 : coordinate system (not supported)
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
        else if (cmd == "PSHELL")
        {
            // The last ditch effort to map PSHELL id to a group name.
            // If ANSA or HMNAME didn't work, it is still possible to
            // have the 'weird' format where the immediately preceeding
            // comment contains the information.

            // Fixed format:
            //  8-16 : pshell id

            label groupId = readLabel(nextNasField(line, linei, 8, freeFormat));

            if (lastComment.size() > 1 && !nameLookup.contains(groupId))
            {
                word groupName = word::validate(lastComment.substr(1));

                if (!groupName.empty())
                {
                    DebugInfo
                        << "PSHELL:" << groupId
                        << " = " << groupName << nl;

                    nameLookup.emplace(groupId, std::move(groupName));
                }
            }
        }
        else if (unhandledCmd.insert(cmd))
        {
            InfoErr
                << "Unhandled Nastran command " << line << nl
                << "File:" << is.name() << " line:" << is.lineNumber()
                << nl;
        }

        // Discard buffered comment (from weird format...)
        lastComment.clear();
    }


    //    Info<< "Read faces:" << dynFaces.size()
    //        << " points:" << dynPoints.size()
    //        << endl;

    if (ignoreElemId)
    {
        dynElemId.clear();
    }

    // Transfer to normal lists
    this->storedPoints().transfer(dynPoints);

    dynFaces.shrink();

    // Build inverse mapping (NASTRAN pointId -> index)
    Map<label> mapPointId(invertToMap(pointId));
    pointId.clearStorage();

    // Relabel faces
    // ~~~~~~~~~~~~~
    for (Face& f : dynFaces)
    {
        for (label& vert : f)
        {
            vert = mapPointId[vert];
        }
    }
    mapPointId.clear();

    DebugInfo
        << "PSHELL names:" << nameLookup << nl;

    // Create default zone names, or from ANSA/Hypermesh information
    List<word> names(dynSizes.size());
    forAllConstIters(zoneLookup, iter)
    {
        const label groupId = iter.key();
        const label zoneId  = iter.val();

        const auto iterName = nameLookup.cfind(groupId);
        if (iterName.good())
        {
            names[zoneId] = iterName.val();
        }
        else
        {
            names[zoneId] = surfZone::defaultName(zoneId);
        }
    }

    this->sortFacesAndStore(dynFaces, dynZones, dynElemId, sorted);

    // Add zones (retaining empty ones)
    this->addZones(dynSizes, names);
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::NASsurfaceFormat<Face>::write
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

    // For simplicity, use fieldFormat::FREE throughout
    fileFormats::NASCore::setPrecision(os, fieldFormat::FREE);

    os  << "CEND" << nl
        << "TITLE = " << os.name().stem() << nl;

    // Print zone names as comment
    forAll(zones, zonei)
    {
        // HYPERMESH extension
        os  << "$HMNAME COMP" << setw(20) << (zonei+1)
            << '"' << zones[zonei].name() << '"' << nl;
    }

    // Write vertex coords with 1-based point Id
    os  << "$ GRID POINTS" << nl
        << "BEGIN BULK" << nl;

    label pointId = 0;
    for (const point& pt : pointLst)
    {
        os  << "GRID" << ','
            << ++pointId << ','
            << 0 << ','  // global coordinate system
            << pt.x() << ',' << pt.y() << ',' << pt.z() << nl;
    }

    os << "$ ELEMENTS" << nl;

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

            elemId = writeShell(os, f, elemId, zoneIndex);
        }

        ++zoneIndex;
    }

    os << "ENDDATA" << nl;
}


// ************************************************************************* //
