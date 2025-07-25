/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::fileFormats::TRIReader

Description
    TRI (triangle) file reader.

    For TRI format (eg, AC3D).
    Each input line has 9 floats (3 points, each 3 floats) followed by hex
    colour. The colour is used to build regions numbered from 0 upwards.

    Reading and stitching similar to the STLReader.

SourceFiles
    TRIReader.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_TRIReader_H
#define Foam_TRIReader_H

#include "STLpoint.H"
#include "labelList.H"
#include "surfaceFormatsCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::TRIReader Declaration
\*---------------------------------------------------------------------------*/

class TRIReader
:
    public surfaceFormatsCore
{
    // Private Data

        bool sorted_;

        //- The points supporting the facets
        List<STLpoint> points_;

        //- The zones associated with the faces
        List<label> zoneIds_;

        //- The zone names, in the order of their first appearance
        List<word> names_;

        //- The solid count, in the order of their first appearance
        List<label> sizes_;


    // Private Member Functions

        bool readFile(const fileName& filename);


public:

    // Constructors

        //- Read from file, filling in the information
        explicit TRIReader(const fileName& filename);


    //- Destructor
    ~TRIReader() = default;


    // Member Functions

        //- Flush all values
        void clear();

        //- Calculate merge points mapping, return old to new pointMap.
        //  Use merge tolerance as per STL ascii
        //  \return number of unique points
        label mergePointsMap(labelList& pointMap) const;

        //- Calculate merge points mapping, return old to new pointMap.
        //  \return number of unique points
        label mergePointsMap(const scalar mergeTol, labelList& pointMap) const;


        //- File read was already in sorted order
        bool is_sorted() const noexcept { return sorted_; }

        //- Return full access to the points
        List<STLpoint>& points() noexcept { return points_; }

        //- Return full access to the zones
        List<label>& zoneIds() noexcept { return zoneIds_; }

        //- The list of solid names in the order of their first appearance
        List<word>& names() noexcept { return names_; }

        //- The list of zone sizes in the order of their first appearance
        List<label>& sizes() noexcept { return sizes_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
