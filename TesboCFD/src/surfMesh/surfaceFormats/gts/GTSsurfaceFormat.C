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

#include "GTSsurfaceFormat.H"
#include "surfaceFormatsCore.H"
#include "clock.H"
#include "Fstream.H"
#include "StringStream.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::GTSsurfaceFormat<Face>::checkIfTriangulated
(
    const UList<Face>& faceLst
)
{
    label nNonTris = 0;

    if (!faceTraits<Face>::isTri())
    {
        for (const auto& f : faceLst)
        {
            if (f.size() != 3)
            {
                ++nNonTris;
            }
        }
    }

    if (nNonTris)
    {
        FatalErrorInFunction
            << "Surface has " << nNonTris << '/' << faceLst.size()
            << " non-triangulated faces - not writing!" << endl;
    }

    return nNonTris == 0;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::GTSsurfaceFormat<Face>::GTSsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::GTSsurfaceFormat<Face>::read
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

    // Read header
    string line = this->getLineNoComment(is);

    label nPoints, nEdges, nElems;
    {
        IStringStream lineStream(line);
        lineStream
            >> nPoints
            >> nEdges
            >> nElems;
    }


    // Write directly into the lists
    auto& pointLst = this->storedPoints();
    auto& faceLst  = this->storedFaces();
    auto& zoneIds  = this->storedZoneIds();

    pointLst.resize(nPoints);
    faceLst.resize(nElems);
    zoneIds.resize(nElems);

    // Read points
    forAll(pointLst, pointi)
    {
        scalar x, y, z;
        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> x >> y >> z;
        }

        pointLst[pointi] = point(x, y, z);
    }

    // Read edges (OpenFOAM indexing)
    edgeList edges(nEdges);
    forAll(edges, edgei)
    {
        label beg, end;
        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> beg >> end;
        }
        edges[edgei] = edge(beg - 1, end - 1);
    }


    // Read triangles. Convert references to edges into pointlabels
    label maxZone = 0;
    forAll(faceLst, facei)
    {
        label e0Label, e1Label, e2Label;
        label zoneI = 0;

        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> e0Label >> e1Label >> e2Label;

            // Optional zone number: read first, then check stream state
            if (lineStream)
            {
                label num;
                lineStream >> num;
                if (!lineStream.bad())
                {
                    zoneI = num;
                    if (maxZone < zoneI)
                    {
                        maxZone = zoneI;
                    }
                }
            }
        }

        // Determine ordering of edges e0, e1
        //  common: common vertex, shared by e0 and e1
        //  e0Far:  vertex on e0 which is not common
        //  e1Far:  vertex on e1 which is not common
        const edge& e0 = edges[e0Label - 1];
        const edge& e1 = edges[e1Label - 1];
        const edge& e2 = edges[e2Label - 1];

        label common01 = e0.commonVertex(e1);
        if (common01 == -1)
        {
            FatalErrorInFunction
                << "Edges 0 and 1 of triangle " << facei
                << " do not share a point.\n"
                << "    edge0:" << e0 << nl
                << "    edge1:" << e1
                << exit(FatalError);
        }

        const label e0Far = e0.otherVertex(common01);
        const label e1Far = e1.otherVertex(common01);

        const label common12 = e1.commonVertex(e2);
        if (common12 == -1)
        {
            FatalErrorInFunction
                << "Edges 1 and 2 of triangle " << facei
                << " do not share a point.\n"
                << "    edge1:" << e1 << nl
                << "    edge2:" << e2
                << exit(FatalError);
        }
        const label e2Far = e2.otherVertex(common12);

        // Does edge2 sit between edge1 and 0?
        if (common12 != e1Far || e2Far != e0Far)
        {
            FatalErrorInFunction
                << "Edges of triangle " << facei
                << " reference more than three points.\n"
                << "    edge0:" << e0 << nl
                << "    edge1:" << e1 << nl
                << "    edge2:" << e2 << nl
                << exit(FatalError);
        }

        faceLst[facei] = Face{e0Far, common01, e1Far};
        zoneIds[facei] = zoneI;
    }


    List<surfZoneIdentifier> newZones(maxZone+1);
    forAll(newZones, zonei)
    {
        newZones[zonei] = surfZoneIdentifier
        (
            surfZoneIdentifier::defaultName(zonei),
            zonei
        );
    }

    this->storedZoneToc().transfer(newZones);
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::GTSsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurface<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary&
)
{
    // ASCII only, allow output compression
    streamOpt.format(IOstreamOption::ASCII);

    const UList<point>& pointLst = surf.points();
    const UList<Face>& faceLst = surf.surfFaces();

    const surfZoneList zones =
    (
        surf.surfZones().size()
      ? surf.surfZones()
      : surfaceFormatsCore::oneZone(faceLst)
    );

    checkIfTriangulated(faceLst);

    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }


    // Write header, print zone names as comment
    os  << "# GTS file" << nl
        << "# Zones:" << nl;

    forAll(zones, zonei)
    {
        os  << "#     " << zonei << "    "
            << zones[zonei].name() << nl;
    }
    os  << "#" << nl;

    os  << "# nPoints  nEdges  nTriangles" << nl
        << pointLst.size() << ' ' << surf.nEdges() << ' '
        << surf.size() << nl;


    // Write vertex coords
    for (const point& pt : pointLst)
    {
        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }


    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = surf.edges();
    const labelList& meshPts = surf.meshPoints();

    for (const edge& e : es)
    {
        os  << meshPts[e.start()] + 1 << ' '
            << meshPts[e.end()] + 1 << nl;
    }

    // Write faces in terms of edges
    const labelListList& faceEs = surf.faceEdges();

    label faceIndex = 0;
    label zoneIndex = 0;

    for (const surfZone& zone : zones)
    {
        for (label nLocal = zone.size(); nLocal--; ++faceIndex)
        {
            const label facei = faceIndex;

            const labelList& fEdges = faceEs[facei];

            os  << fEdges[0] + 1 << ' '
                << fEdges[1] + 1 << ' '
                << fEdges[2] + 1 << ' '
                << zoneIndex << nl;
        }

        ++zoneIndex;
    }
}


template<class Face>
void Foam::fileFormats::GTSsurfaceFormat<Face>::write
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
    const UList<Face>& faceLst = surf.surfFaces();
    const UList<label>& zoneIds = surf.zoneIds();
    const UList<surfZoneIdentifier>& zoneToc = surf.zoneToc();

    checkIfTriangulated(faceLst);

    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }


    // Write header, print zone names as comment
    os  << "# GTS file" << nl
        << "# Zones:" << nl;

    forAll(zoneToc, zonei)
    {
        os  << "#     " << zonei << "    "
            << zoneToc[zonei].name() << nl;
    }
    os  << "#" << nl;

    os  << "# nPoints  nEdges  nTriangles" << nl
        << pointLst.size() << ' ' << surf.nEdges() << ' '
        << surf.size() << nl;


    // Write vertex coords
    for (const point& pt : pointLst)
    {
        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }


    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = surf.edges();
    const labelList& meshPts = surf.meshPoints();

    for (const edge& e : es)
    {
        os  << meshPts[e.start()] + 1 << ' '
            << meshPts[e.end()] + 1 << nl;
    }

    // Write faces in terms of edges.
    const labelListList& faceEs = surf.faceEdges();

    forAll(faceLst, facei)
    {
        const labelList& fEdges = faceEs[facei];

        os  << fEdges[0] + 1 << ' '
            << fEdges[1] + 1 << ' '
            << fEdges[2] + 1 << ' '
            << zoneIds[facei] << nl;
    }
}


// ************************************************************************* //
