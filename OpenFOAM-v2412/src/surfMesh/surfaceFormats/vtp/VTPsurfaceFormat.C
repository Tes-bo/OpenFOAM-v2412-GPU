/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "VTPsurfaceFormat.H"
#include <fstream>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::writePolys
(
    vtk::formatter& format,
    const UList<Face>& faces
)
{
    format.tag(vtk::fileTag::POLYS);

    //
    // 'connectivity'
    //
    {
        label nVerts = 0;
        for (const auto& f : faces)
        {
            nVerts += f.size();
        }

        const uint64_t payLoad = vtk::sizeofData<label>(nVerts);

        format.beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
        format.writeSize(payLoad);

        for (const auto& f : faces)
        {
            vtk::writeList(format, f);
        }

        format.flush();
        format.endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        const uint64_t payLoad = vtk::sizeofData<label>(faces.size());

        format.beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
        format.writeSize(payLoad);

        label off = 0;
        for (const auto& f : faces)
        {
            off += f.size();

            format.write(off);
        }

        format.flush();
        format.endDataArray();
    }

    format.endTag(vtk::fileTag::POLYS);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    IOstreamOption,
    const dictionary& options
)
{
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

    vtk::outputOptions opts = formatOptions(options);

    std::ofstream os(filename, std::ios::binary);

    autoPtr<vtk::formatter> format = opts.newFormatter(os);

    writeHeader(format(), pointLst, faceLst.size());

    if (useFaceMap)
    {
        format().tag(vtk::fileTag::POLYS);

        //
        // 'connectivity'
        //
        {
            label nVerts = 0;
            for (const auto& f : faceLst)
            {
                nVerts += f.size();
            }

            const uint64_t payLoad = vtk::sizeofData<label>(nVerts);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            label faceIndex = 0;
            for (const surfZone& zone : zones)
            {
                forAll(zone, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];

                    vtk::writeList(format(), f);
                }
            }

            format().flush();
            format().endDataArray();
        }


        //
        // 'offsets'  (connectivity offsets)
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(faceLst.size());

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            label off = 0, faceIndex = 0;
            for (const surfZone& zone : zones)
            {
                forAll(zone, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];

                    off += f.size();

                    format().write(off);
                }
            }

            format().flush();
            format().endDataArray();
        }

        format().endTag(vtk::fileTag::POLYS);
    }
    else
    {
        // Easy to write polys without a faceMap
        writePolys(format(), faceLst);
    }

    // Write regions (zones) as CellData
    if (zones.size() > 1)
    {
        writeCellData(format(), zones);
    }

    writeFooter(format());
}


template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    IOstreamOption,
    const dictionary& options
)
{
    vtk::outputOptions opts = formatOptions(options);

    std::ofstream os(filename, std::ios::binary);

    autoPtr<vtk::formatter> format = opts.newFormatter(os);

    const UList<Face>& faceLst = surf.surfFaces();

    writeHeader(format(), surf.points(), faceLst.size());

    // Easy to write polys without a faceMap
    writePolys(format(), faceLst);

    // Write regions (zones) as CellData
    writeCellData(format(), surf.zoneIds());

    writeFooter(format());
}


// ************************************************************************* //
