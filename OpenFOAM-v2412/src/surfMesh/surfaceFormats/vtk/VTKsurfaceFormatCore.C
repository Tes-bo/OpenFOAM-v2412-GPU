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

#include "VTKsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vtk::outputOptions
Foam::fileFormats::VTKsurfaceFormatCore::formatOptions
(
    const dictionary& dict,
    vtk::outputOptions opts
)
{
    opts.legacy(true);  // Legacy. Use VTPsurfaceFormat for non-legacy
    opts.append(false); // No append format for legacy

    opts.ascii
    (
        IOstreamOption::ASCII
     == IOstreamOption::formatEnum("format", dict, IOstreamOption::ASCII)
    );

    opts.precision
    (
        dict.getOrDefault("precision", IOstream::defaultPrecision())
    );

    return opts;
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeHeader
(
    vtk::formatter& format,
    const UList<point>& pts
)
{
    vtk::legacy::fileHeader<vtk::fileTag::POLY_DATA>
    (
        format,
        ("surface written " + clock::dateTime())
    );

    vtk::legacy::beginPoints(format.os(), pts.size());

    vtk::writeList(format, pts);
    format.flush();
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeCellData
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

    vtk::legacy::beginCellData(format, nFaces, 1);      // 1 field
    vtk::legacy::intField<1>(format, "region", nFaces); // 1 component

    label zoneId = 0;
    for (const surfZone& z : zones)
    {
        vtk::write(format, zoneId, z.size());
        ++zoneId;
    }
    format.flush();
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeCellData
(
    vtk::formatter& format,
    const labelUList& zoneIds
)
{
    // Zone ids as CellData

    // Number of faces
    const label nFaces = zoneIds.size();

    vtk::legacy::beginCellData(format, nFaces, 1);      // 1 field
    vtk::legacy::intField<1>(format, "region", nFaces); // 1 component

    vtk::writeList(format, zoneIds);
    format.flush();
}


// ************************************************************************* //
