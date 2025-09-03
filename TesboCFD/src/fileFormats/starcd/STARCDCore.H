/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::fileFormats::STARCDCore

Description
    Core routines used when reading/writing PROSTAR vrt/cel/bnd files.

SourceFiles
    STARCDCore.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDCore_H
#define STARCDCore_H

#include "IFstream.H"
#include "Enum.H"
#include "Map.H"
#include "point.H"
#include "FixedList.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::STARCDCore Declaration
\*---------------------------------------------------------------------------*/

class STARCDCore
{
public:

    // Public Data, Declarations

        //- Enumeration defining the file headers
        enum fileHeader
        {
            HEADER_CEL,
            HEADER_VRT,
            HEADER_BND
        };

        //- Enumeration defining the file extensions
        enum fileExt
        {
            CEL_FILE,
            VRT_FILE,
            BND_FILE,
            INP_FILE
        };

        //- Basic material type for STARCD/PROSTAR files
        enum matlType
        {
            starcdFluidType   = 1,
            starcdSolidType   = 2,
            starcdBaffleType  = 3,
            starcdShellType   = 4,
            starcdLineType    = 5,
            starcdPointType   = 6
        };

        //- Shape-Type for STARCD/PROSTAR files
        enum shapeType
        {
            starcdPoint = 1,
            starcdLine  = 2,
            starcdShell = 3,
            starcdHex   = 11,
            starcdPrism = 12,
            starcdTet   = 13,
            starcdPyr   = 14,
            starcdPoly  = 255
        };


        //- The name for default (unassigned) boundaries
        static const char* const defaultBoundaryName;

        //- The name we have chosen for default (unassigned) solid boundaries.
        //  Slightly distinguished from the regular default name.
        static const char* const defaultSolidBoundaryName;


private:

    // Private Data

        static const Enum<fileHeader> fileHeaders_;
        static const Enum<fileExt>    fileExtensions_;


protected:

    // Protected Member Functions

        //- Face addressing from PROSTAR faces to OpenFOAM faces.
        //  For hex, prism, tet, pyr primitive shapes.
        static const Map<FixedList<int, 6>> starToFoamFaceAddr;

        //- Face addressing from OpenFOAM faces to PROSTAR faces.
        //  For hex, prism, tet, pyr primitive shapes.
        static const Map<FixedList<int, 6>> foamToStarFaceAddr;


    // Constructors

        //- Construct null
        STARCDCore() = default;

        //- Read header and check signature PROSTAR_(CELL|VERTEX|BOUNDARY)
        static bool readHeader(IFstream& is, const enum fileHeader header);

        //- Write header for fileType (CELL|VERTEX|BOUNDARY)
        static void writeHeader(Ostream& os, const enum fileHeader header);


public:

    // Public Member Functions

        //- Resolve base file-name for the given file-type
        static fileName starFileName
        (
            const fileName& baseName,
            const enum fileExt ext
        );


        //- Remove existing PROSTAR files for the given base file-name
        static void removeFiles(const fileName& baseName);


        //- Read points from a (.vrt) file, return the max prostar id used.
        //
        // The file format is as follows:
        // \verbatim
        // Line 1:
        //   PROSTAR_VERTEX  newline
        //
        // Line 2:
        //   {version} 0 0 0 0 0 0 0  newline
        //
        // Body:
        //   {vertexId}  {x}  {y}  {z}  newline
        // \endverbatim
        static label readPoints
        (
            IFstream& is,
            List<point>& points,
            List<label>& ids
        );

        //- Write header and points to (.vrt) file, optionally with scaling
        static void writePoints
        (
            Ostream& os,
            const UList<point>& points,
            const scalar scaleFactor = 1.0
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
