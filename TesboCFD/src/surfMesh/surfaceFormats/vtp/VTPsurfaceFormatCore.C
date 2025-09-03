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

#include "VTPsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vtk::outputOptions
Foam::fileFormats::VTPsurfaceFormatCore::formatOptions
(
    const dictionary& dict,
    vtk::outputOptions opts
)
{
    opts.legacy(false); // Non-legacy. Use VTKsurfaceFormat for legacy
    opts.append(false); // No append format

    opts.ascii
    (
        IOstreamOption::ASCII
     == IOstreamOption::formatEnum("format", dict, IOstreamOption::BINARY)
    );

    opts.precision
    (
        dict.getOrDefault("precision", IOstream::defaultPrecision())
    );

    return opts;
}


void Foam::fileFormats::VTPsurfaceFormatCore::writeHeader
(
    vtk::formatter& format,
    const UList<point>& pts,
    const label nFaces
)
{
    // XML (inline)

    format
        .xmlHeader()
        .xmlComment("surface written " + clock::dateTime())
        .beginVTKFile<vtk::fileTag::POLY_DATA>();

    // <Piece>
    format
        .tag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, pts.size(),
            vtk::fileAttr::NUMBER_OF_POLYS,  nFaces
        );


    // Points

    const uint64_t payLoad = vtk::sizeofData<float, 3>(pts.size());

    format.tag(vtk::fileTag::POINTS)
        .beginDataArray<float, 3>(vtk::dataArrayAttr::POINTS);

    format.writeSize(payLoad);

    vtk::writeList(format, pts);
    format.flush();

    format
        .endDataArray()
        .endTag(vtk::fileTag::POINTS);
}


void Foam::fileFormats::VTPsurfaceFormatCore::writeFooter
(
    vtk::formatter& format
)
{
    format.endPiece();  //<-- slight cheat. </Piece> too

    format.endTag(vtk::fileTag::POLY_DATA)
        .endVTKFile();
}


void Foam::fileFormats::VTPsurfaceFormatCore::writeCellData
(
    vtk::formatter& format,
    const UList<surfZone>& zones
)
{
    // Zone ids as CellData

    // Number of faces covered by the zones
    label nFaces = 0;
    for (const surfZone& z : zones)
    {
        nFaces += z.size();
    }

    const uint64_t payLoad = vtk::sizeofData<label>(nFaces);

    format.beginCellData();
    format.beginDataArray<label>("region");
    format.writeSize(payLoad);

    label zoneId = 0;
    for (const surfZone& z : zones)
    {
        vtk::write(format, zoneId, z.size());
        ++zoneId;
    }

    format.flush();
    format.endDataArray();

    format.endCellData();
}


void Foam::fileFormats::VTPsurfaceFormatCore::writeCellData
(
    vtk::formatter& format,
    const labelUList& zoneIds
)
{
    // Zone ids as CellData

    const uint64_t payLoad = vtk::sizeofData<label>(zoneIds.size());

    format.beginCellData();
    format.beginDataArray<label>("region");

    format.writeSize(payLoad);
    vtk::writeList(format, zoneIds);

    format.flush();
    format.endDataArray();

    format.endCellData();
}


// ************************************************************************* //
