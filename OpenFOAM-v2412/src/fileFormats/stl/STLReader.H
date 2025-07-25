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
    Foam::fileFormats::STLReader

Description
    Internal class used by the STLsurfaceFormat and triSurface.

SourceFiles
    STLReader.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef Foam_STLReader_H
#define Foam_STLReader_H

#include "STLCore.H"
#include "labelledTri.H"
#include "IFstream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Detail
{

// Forward declarations
class STLAsciiParse;

} // End namespace Detail


namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::STLReader Declaration
\*---------------------------------------------------------------------------*/

class STLReader
:
    public STLCore
{
    // Private Data

        bool sorted_;

        //- The points supporting the facets
        List<STLpoint> points_;

        //- The zones associated with the faces
        List<label> zoneIds_;

        //- The solid names, in the order of their first appearance
        List<word> names_;

        //- The solid count, in the order of their first appearance
        List<label> sizes_;

        //- The STL format used
        STLFormat format_;


    // Private Member Functions

        //- Transfer parsed information to normal lists
        void transfer(Detail::STLAsciiParse& parsed);

        //- Parse/read ASCII using Flex-based parser
        bool readAsciiFlex(const fileName& filename);

        //- Parse/read ASCII using Ragel-based parser
        bool readAsciiRagel(const fileName& filename);

        //- Parse/read ASCII using simple handwritten parser
        bool readAsciiManual(const fileName& filename);

        //- Parse/read ASCII
        bool readASCII(const fileName& filename);

        //- Read BINARY
        bool readBINARY(const fileName& filename);

        //- Read ASCII or BINARY
        bool readFile(const fileName& filename, const STLFormat format);


        //- No copy construct
        STLReader(const STLReader&) = delete;

        //- No copy assignment
        void operator=(const STLReader&) = delete;


public:

    // Static Data

        //- ASCII parser types (0=Flex, 1=Ragel, 2=Manual)
        static int parserType;


    // Constructors

        //- Read from file, filling in the information.
        //  Auto-detect ASCII/BINARY format.
        explicit STLReader(const fileName& filename);

        //- Read from file, filling in the information.
        //  Manually selected choice of ASCII/BINARY/UNKNOWN(detect) formats.
        STLReader(const fileName& filename, const STLFormat format);


    //- Destructor
    ~STLReader() = default;


    // Member Functions

        //- Flush all values
        void clear();

        //- Calculate merge points mapping, return old to new pointMap.
        //  The merge tolerance based on ASCII or BINARY input format.
        //  \return number of unique points
        label mergePointsMap(labelList& pointMap) const;

        //- Calculate merge points mapping, return old to new pointMap.
        //  \return number of unique points
        label mergePointsMap(const scalar mergeTol, labelList& pointMap) const;

        //- File read was already sorted?
        bool is_sorted() const noexcept { return sorted_; }

        //- Return full access to the points
        List<STLpoint>& points() noexcept { return points_; }

        //- Return full access to the zoneIds
        List<label>& zoneIds() noexcept { return zoneIds_; }

        //- The list of solid names in the order of their first appearance
        List<word>& names() noexcept { return names_; }

        //- The list of solid sizes in the order of their first appearance
        List<label>& sizes() noexcept { return sizes_; }

        //- The STL format used (ASCII or BINARY)
        enum STLFormat stlFormat() const noexcept { return format_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
