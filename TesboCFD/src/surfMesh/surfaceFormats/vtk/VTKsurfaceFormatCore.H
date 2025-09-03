/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

Class
    Foam::fileFormats::VTKsurfaceFormatCore

Description
    Internal class used by the VTKsurfaceFormat
    Format is LEGACY_ASCII

    \heading Output Options
    \table
        Property    | Description                       | Required | Default
        format      | ascii or binary format            | no | ascii
        precision   | Write precision in ascii          | no | same as IOstream
    \endtable

SourceFiles
    VTKsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef VTKsurfaceFormatCore_H
#define VTKsurfaceFormatCore_H

#include "point.H"
#include "surfZone.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
              Class fileFormats::VTKsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class VTKsurfaceFormatCore
{
protected:

    // Protected Static Member Functions

        //- Extract format options (default format LEGACY_ASCII)
        static vtk::outputOptions formatOptions
        (
            const dictionary& dict,
            vtk::outputOptions opts = vtk::formatType::LEGACY_ASCII
        );

        //- Write header information with points
        static void writeHeader
        (
            vtk::formatter& format,
            const UList<point>& pts
        );

        //- Write regions (zones) information as CellData
        static void writeCellData
        (
            vtk::formatter& format,
            const UList<surfZone>& zones
        );

        //- Write regions (zones) information as CellData
        static void writeCellData
        (
            vtk::formatter& format,
            const labelUList& zoneIds
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
