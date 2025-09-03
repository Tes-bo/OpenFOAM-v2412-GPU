/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::fileFormats::STLCore

Description
    Core routines used when reading/writing STL files.

SourceFiles
    STLCore.C

\*---------------------------------------------------------------------------*/

#ifndef STLCore_H
#define STLCore_H

#include "STLpoint.H"
#include "STLtriangle.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class fileFormats::STLCore Declaration
\*---------------------------------------------------------------------------*/

class STLCore
{
public:

    // Public Data Types

        //- Enumeration for the format of data in the stream
        enum STLFormat
        {
            ASCII,  //!< ASCII
            BINARY, //!< BINARY
            UNKNOWN //!< Detect based on (input) content or (output) extension
        };


protected:

    // Protected Member Functions

        //- Detect 'stlb' extension as binary when format = UNKNOWN.
        //  Otherwise test if format == BINARY.
        static bool isBinaryName
        (
            const fileName& filename,
            const STLFormat format
        );


        //- Check contents to detect if the file is a binary STL.
        //  Return the estimated number of triangles or 0 on error.
        static int detectBinaryHeader(const fileName& filename);


        //- Read STL binary file header.
        //  Return the opened file stream and estimated number of triangles.
        //  The stream is invalid and number of triangles is 0 on error.
        static std::unique_ptr<std::istream>
        readBinaryHeader
        (
            const fileName& filename,
            label& nTrisEstimated
        );

        //- Write STL binary file and number of triangles to stream
        static void writeBinaryHeader(ostream& os, uint32_t nTris);


    // Constructors

        //- Default construct
        STLCore() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
