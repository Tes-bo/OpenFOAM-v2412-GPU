/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::cell

Description
    A cell is defined as a list of faces with extra functionality

SourceFiles
    cellI.H
    cell.C
    oppositeCellFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cell_H
#define Foam_cell_H

#include "faceList.H"
#include "oppositeFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class cell Declaration
\*---------------------------------------------------------------------------*/

class cell
:
    public labelList
{
public:

    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Default construct
        constexpr cell() noexcept = default;

        //- Construct given size, with invalid point labels (-1)
        inline explicit cell(const label sz);

        //- Copy construct from list of face labels
        inline explicit cell(const labelUList& list);

        //- Move construct from list of face labels
        inline explicit cell(labelList&& list);

        //- Copy construct from list of face labels
        template<unsigned N>
        inline explicit cell(const FixedList<label, N>& list);

        //- Copy construct from subset of face labels
        inline cell(const labelUList& list, const labelUList& indices);

        //- Copy construct from subset of face labels
        template<unsigned N>
        inline cell(const labelUList& list, const FixedList<label, N>& indices);

        //- Construct from Istream
        inline explicit cell(Istream& is);


    // Member Functions

        //- Return number of faces
        inline label nFaces() const noexcept;

        //- Return unordered list of cell vertices given the list of faces
        labelList labels(const faceUList& meshFaces) const;

        //- Return the cell vertices given the list of faces and mesh points
        pointField points
        (
            const faceUList& meshFaces,
            const UList<point>& meshPoints
        ) const;

        //- Return cell edges
        edgeList edges(const faceUList& meshFaces) const;

        //- Return index of opposite face
        label opposingFaceLabel
        (
            const label masterFaceLabel,
            const faceUList& meshFaces
        ) const;

        //- Return opposite face oriented the same way as the master face
        oppositeFace opposingFace
        (
            const label masterFaceLabel,
            const faceUList& meshFaces
        ) const;


        // Instructions for centre and mag:
        // When one wants to access the cell centre and magnitude, the
        // functionality on the mesh level should be used in preference to the
        // functions provided here. They do not rely to the functionality
        // implemented here, provide additional checking and are more efficient.
        // The cell::centre and cell::mag functions may be removed in the
        // future.

        //- Returns cell centre
        point centre
        (
            const UList<point>& meshPoints,
            const faceUList& meshFaces
        ) const;

        //- Returns cell volume
        scalar mag
        (
            const UList<point>& meshPoints,
            const faceUList& meshFaces
        ) const;

        //- The bounding box for the cell
        Pair<point> box
        (
            const UList<point>& meshPoints,
            const faceUList& meshFaces
        ) const;

        //- The bounding box for the cell associated with given mesh
        Pair<point> box(const primitiveMesh& mesh) const;
};


// Global Operators

//- Test if both cells are the same size and contain the same points
//  The internal point ordering is ignored
bool operator==(const cell& a, const cell& b);

//- Test if the cells differ (different size or different points)
inline bool operator!=(const cell& a, const cell& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
