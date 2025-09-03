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
    Foam::fileFormats::VTPsurfaceFormatCore

Description
    Internal class used by the VTPsurfaceFormat.
    The default format is INLINE_BASE64

    \heading Output Options
    \table
        Property    | Description                       | Required | Default
        format      | ascii or binary format            | no | binary
        precision   | Write precision in ascii          | no | same as IOstream
    \endtable

SourceFiles
    VTPsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef VTPsurfaceFormatCore_H
#define VTPsurfaceFormatCore_H

#include "point.H"
#include "surfZone.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
              Class fileFormats::VTPsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class VTPsurfaceFormatCore
{
protected:

    // Protected Static Member Functions

        //- Extract format options (default format INLINE_BASE64)
        static vtk::outputOptions formatOptions
        (
            const dictionary& dict,
            vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Write header information with points
        static void writeHeader
        (
            vtk::formatter& format,
            const UList<point>& pts,
            const label nFaces
        );

        //- Write footer
        static void writeFooter(vtk::formatter& format);

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
