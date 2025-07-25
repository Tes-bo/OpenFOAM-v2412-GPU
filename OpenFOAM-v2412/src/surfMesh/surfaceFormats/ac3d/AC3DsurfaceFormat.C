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

#include "AC3DsurfaceFormat.H"
#include "StringStream.H"
#include "primitivePatch.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::AC3DsurfaceFormat<Face>::AC3DsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::AC3DsurfaceFormat<Face>::read
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

    string line, cmd, args;

    is.getLine(line);

    // Verify version
    {
        const string version = line.substr(4);

        if (version != "b")
        {
            WarningInFunction
                << "When reading AC3D file " << filename
                << " read header " << line << " with version "
                << version << endl
                << "Only tested reading with version 'b'."
                << " This might give problems" << endl;
        }
    }


    if (!cueTo(is, "OBJECT", args) || args != "world")
    {
        FatalErrorInFunction
            << "Cannot find 'OBJECT world' in file " << filename << nl
            << exit(FatalError);
    }

    // Number of kids is the number of zones
    args = cueToOrDie(is, "kids");
    const label nZones = parse<int>(args);

    // Start of vertices for object/zones
    label vertexOffset = 0;

    DynamicList<point> dynPoints;
    DynamicList<Face>  dynFaces;
    List<word>         names(nZones);
    List<label>        sizes(nZones, Zero);

    for (label zoneI = 0; zoneI < nZones; ++zoneI)
    {
        names[zoneI] = surfZone::defaultName(zoneI);

        args = cueToOrDie(is, "OBJECT", "while reading " + names[zoneI]);

        // number of vertices for this zone
        label  nZonePoints = 0;
        vector location(Zero);
        // tensor rotation(I);

        // Read all info for current zone
        while (is.good())
        {
            // Read line and get first word. If end of file break since
            // zone should always end with 'kids' command ?not sure.
            if (!readCmd(is, cmd, args))
            {
                FatalErrorInFunction
                    << "Did not read up to 'kids 0' while reading zone "
                    << zoneI << " from file " << filename << nl
                    << exit(FatalError);
            }

            if (cmd == "name")
            {
                // name %s
                const string str = parse<string>(args);
                names[zoneI] = word::validate(str);
            }
            else if (cmd == "rot")
            {
                // rot  %f %f %f  %f %f %f  %f %f %f

                // IStringStream lineStream(args);
                //
                // lineStream
                //     >> rotation.xx() >> rotation.xy() >> rotation.xz()
                //     >> rotation.yx() >> rotation.yy() >> rotation.yz()
                //     >> rotation.zx() >> rotation.zy() >> rotation.zz();

                WarningInFunction
                    << "rot (rotation tensor) command not implemented"
                    << "Line:" << cmd << ' ' << args << endl
                    << "while reading zone " << zoneI << endl;
            }
            else if (cmd == "loc")
            {
                // loc  %f %f %f
                IStringStream lineStream(args);

                lineStream
                    >> location.x()
                    >> location.y()
                    >> location.z();
            }
            else if (cmd == "numvert")
            {
                // numvert  %d
                nZonePoints = parse<int>(args);

                for (label vertI = 0; vertI < nZonePoints; ++vertI)
                {
                    is.getLine(line);
                    IStringStream lineStream(line);

                    point pt;
                    lineStream
                        >> pt.x() >> pt.y() >> pt.z();

                    // Offset with current translation vector
                    dynPoints.append(location + pt);
                }
            }
            else if (cmd == "numsurf")
            {
                const label nFaces = parse<int>(args);

                for (label facei = 0; facei < nFaces; ++facei)
                {
                    const string errorMsg =
                        string(" while reading face ")
                            + Foam::name(facei) + " on zone "
                            + Foam::name(zoneI)
                            + " from file " + filename;

                    cueToOrDie(is, "SURF", errorMsg);
                    cueToOrDie(is, "mat", errorMsg);
                    args = cueToOrDie(is, "refs", errorMsg);

                    const label nVert = parse<int>(args);

                    List<label> verts(nVert);
                    forAll(verts, vertI)
                    {
                        is.getLine(line);
                        verts[vertI] = vertexOffset + parse<int>(line);
                    }

                    const labelUList& f = static_cast<const labelUList&>(verts);

                    if (faceTraits<Face>::isTri() && f.size() > 3)
                    {
                        // simple face triangulation about f[0]
                        // points may be incomplete
                        for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
                        {
                            label fp2 = f.fcIndex(fp1);

                            dynFaces.append(Face{f[0], f[fp1], f[fp2]});
                            sizes[zoneI]++;
                        }
                    }
                    else
                    {
                        dynFaces.append(Face(f));
                        sizes[zoneI]++;
                    }
                }

                // Done the current zone.
                // Increment the offset vertices are stored at
                vertexOffset += nZonePoints;
            }
            else if (cmd == "kids")
            {
                // 'kids' denotes the end of the current zone.
                const label nKids = parse<int>(args);

                if (nKids != 0)
                {
                    FatalErrorInFunction
                        << "Can only read objects without kids."
                        << " Encountered " << nKids << " kids when"
                        << " reading zone " << zoneI
                        << exit(FatalError);
                }

                // Done reading current zone
                break;
            }
        }
    }

    // transfer to normal lists
    this->storedPoints().transfer(dynPoints);
    this->storedFaces().transfer(dynFaces);

    // Add zones (retaining empty ones)
    this->addZones(sizes, names);
    this->addZonesToFaces(); // for labelledTri
    this->stitchFaces(SMALL);

    return true;
}


namespace Foam
{
// file-scope writing of a patch of faces
template<class Patch>
static void writeZone
(
    Ostream& os,
    const Patch& patch,
    const word& name,
    const label zoneI
)
{
    // An isolated surface region (patch).
    os  << "OBJECT poly" << nl
        << "name \"" << name << "\"" << nl;

    os << "numvert " << patch.nPoints() << nl;

    for (const point& pt : patch.localPoints())
    {
        os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }

    os << "numsurf " << patch.size() << nl;

    for (const auto& f : patch.localFaces())
    {
        os  << "SURF 0x20" << nl          // polygon
            << "mat " << zoneI << nl
            << "refs " << f.size() << nl;

        for (const label verti : f)
        {
            os << verti << " 0 0" << nl;
        }
    }

    os << "kids 0" << endl;
}
}

template<class Face>
void Foam::fileFormats::AC3DsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary&
)
{
    // ASCII only, allow output compression
    streamOpt.format(IOstreamOption::ASCII);

    const pointField& pointLst = surf.points();
    const UList<Face>& faceLst = surf.surfFaces();

    const surfZoneList zones =
    (
        surf.surfZones().size()
      ? surf.surfZones()
      : surfaceFormatsCore::oneZone(faceLst)
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    writeHeader(os, zones);

    if (zones.size() == 1)
    {
        PrimitivePatch<UList<Face>, const pointField&> patch
        (
            faceLst, pointLst
        );

        writeZone(os, patch, zones[0].name(), 0);
        return;
    }

    label zoneIndex = 0;
    for (const surfZone& zone : zones)
    {
        if (useFaceMap)
        {
            typedef UIndirectList<Face> FaceListType;

            SubList<label> zoneMap(surf.faceMap(), zone.range());

            PrimitivePatch<FaceListType, const pointField&> patch
            (
                FaceListType(faceLst, zoneMap),
                pointLst
            );

            writeZone(os, patch, zone.name(), zoneIndex);
        }
        else
        {
            typedef SubList<Face> FaceListType;

            PrimitivePatch<FaceListType, const pointField&> patch
            (
                FaceListType(faceLst, zone.range()),
                pointLst
            );

            writeZone(os, patch, zone.name(), zoneIndex);
        }

        ++zoneIndex;
    }
}


template<class Face>
void Foam::fileFormats::AC3DsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary&
)
{
    // ASCII only, allow output compression
    streamOpt.format(IOstreamOption::ASCII);

    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    labelList faceMap;
    List<surfZone> zoneLst = surf.sortedZones(faceMap);

    if (zoneLst.size() <= 1)
    {
        const surfZoneList zones =
        (
            zoneLst.size()
          ? zoneLst
          : surfaceFormatsCore::oneZone(surf.surfFaces())
        );

        writeHeader(os, zones);
        writeZone(os, surf, zones[0].name(), 0);
        return;
    }

    writeHeader(os, zoneLst);

    label zoneIndex = 0;
    for (const surfZone& zone : zoneLst)
    {
        typedef UIndirectList<Face> FaceListType;

        SubList<label> zoneMap(faceMap, zone.range());

        PrimitivePatch<FaceListType, const pointField&> patch
        (
            FaceListType(surf.surfFaces(), zoneMap),
            surf.points()
        );

        writeZone(os, patch, zone.name(), zoneIndex);

        ++zoneIndex;
    }
}


// ************************************************************************* //
