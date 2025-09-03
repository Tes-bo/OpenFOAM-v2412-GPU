/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::fileFormats::NASCore

Description
    Core routines used when reading/writing NASTRAN files.

SourceFiles
    NASCore.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileFormats_NASCore_H
#define Foam_fileFormats_NASCore_H

#include "scalar.H"
#include "Enum.H"
#include "face.H"
#include "point.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class fileFormats::NASCore Declaration
\*---------------------------------------------------------------------------*/

class NASCore
{
public:

    // Public Data

        //- Named/registered debug switch: 'nastran.debug'
        static int debug;


        //- File field formats
        enum fieldFormat
        {
            SHORT,      //!< Short format (field width = 8)
            LONG,       //!< Long format (field width = 16)
            FREE        //!< Free format (comma-separated fields)
        };

        //- Selection names for the NASTRAN file field formats
        static const Enum<fieldFormat> fieldFormatNames;

        //- Output load format
        enum loadFormat
        {
            PLOAD2,     //!< Face load (eg, pressure)
            PLOAD4      //!< Vertex load
        };

        //- Selection names for the NASTRAN load formats
        static const Enum<loadFormat> loadFormatNames;


    // Constructors

        //- Default construct
        NASCore() noexcept = default;


    // Public Static Member Functions

        //- Extract numbers from things like "-2.358-8" (same as "-2.358e-8")
        static scalar readNasScalar(const std::string& str);

        //- A std::string::substr() variant to handle fixed-format and
        //- free-format NASTRAN.
        //  Returns the substr until the next comma (if found)
        //  or the given fixed width
        static std::string nextNasField
        (
            //! The string to extract from
            const std::string& str,
            //! [in,out] The parse position within \p str
            std::string::size_type& pos,
            //! The fixed-format width to use (if comma is not found)
            const std::string::size_type width,
            //! The input is known to be free-format
            const bool free_format = false
        );


        //- Set output stream precision and format flags
        static void setPrecision(Ostream& os, const fieldFormat format);


        //- Write initial keyword (eg, 'GRID' or 'GRID*') followed by the
        //- requisite number of spaces for the field-width
        static Ostream& writeKeyword
        (
            Ostream& os,
            const word& keyword,
            const fieldFormat format
        );

        //- Write a GRID point
        static void writeCoord
        (
            Ostream& os,
            const point& p,
            const label pointId,        //!< 0-based Point Id
            const fieldFormat format
        );

        //- Calculate face decomposition for non tri/quad faces
        //
        //  \param points the surface points
        //  \param faces  the surface faces
        //  \param decompOffsets begin/end offsets (size+1) into decompFaces
        //  \param decompFaces  List of non-tri/quad decomposed into triangles
        //
        //  \return number of decomposed faces
        static label faceDecomposition
        (
            const UList<point>& points,
            const UList<face>& faces,
            labelList& decompOffsets,
            DynamicList<face>& decompFaces
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
