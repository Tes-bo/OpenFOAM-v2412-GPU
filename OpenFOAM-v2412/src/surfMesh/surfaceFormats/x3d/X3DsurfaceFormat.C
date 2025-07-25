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

#include "X3DsurfaceFormat.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::X3DsurfaceFormat<Face>::write
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

    // for no zones, suppress the group name
    const surfZoneList zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst, word::null)
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

    writeHeader(os);
    beginGroup(os);
    writeAppearance(os);

    // NOTE: we could provide an optimized IndexedTriangleSet output for
    // triangulated surfaces too

    os  <<
        "  <IndexedFaceSet coordIndex='\n";

    label faceIndex = 0;
    for (const surfZone& zone : zones)
    {
        const label nLocalFaces = zone.size();

        if (useFaceMap)
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                for (const label vrti : f)
                {
                    os << vrti << ' ';
                }
                os << "-1\n";
            }
        }
        else
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceIndex++];

                for (const label vrti : f)
                {
                    os << vrti << ' ';
                }
                os << "-1\n";
            }
        }
    }

    os  <<
        "' >\n";

    writePoints(os, pointLst);

    os  <<
        "   </IndexedFaceSet>\n";

    endGroup(os);
    writeFooter(os);
}


// ************************************************************************* //
