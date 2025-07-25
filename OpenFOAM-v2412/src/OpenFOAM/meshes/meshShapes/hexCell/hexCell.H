/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::hexCell

Description
    A hexahedral cell primitive.

    It is important that the face/edge ordering is identical for
    a hexahedral cell shape model and a hexCell

SourceFiles
    hexCell.C
    hexCellI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_hexCell_H
#define Foam_hexCell_H

#include "FixedList.H"
#include "faceList.H"
#include "edgeList.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellShape;

/*---------------------------------------------------------------------------*\
                           Class hexCell Declaration
\*---------------------------------------------------------------------------*/

class hexCell
:
    public FixedList<label, 8>
{
    // Static Data Members

        //- The model faces for HEX
        static const label modelFaces_[6][4];

        //- The model edges for HEX
        static const label modelEdges_[12][2];


    // Private Member Functions

        //- Copy vertices for given face - no checks
        inline void copyFace(Foam::face& f, const label facei) const;


public:

    // Generated Methods

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Default construct, with invalid point labels (-1)
        inline hexCell();

        //- Construct from an initializer list of eight point labels
        inline explicit hexCell(std::initializer_list<label> list);

        //- Construct from FixedList of eight point labels
        inline hexCell(const FixedList<label, 8>& list);

        //- Copy construct from a subset of point labels
        inline hexCell
        (
            const labelUList& list,
            const FixedList<label, 8>& indices
        );

        //- Copy construct from a subset of point labels
        template<unsigned AnyNum>
        inline hexCell
        (
            const FixedList<label, AnyNum>& list,
            const FixedList<label, 8>& indices
        );

        //- Construct from Istream
        inline explicit hexCell(Istream& is);


    // Member Functions

        //- Number of points for HEX
        static constexpr label nPoints() noexcept
        {
            return 8;
        }

        //- Number of edges for HEX
        static constexpr label nEdges() noexcept
        {
            return 12;
        }

        //- Number of faces for HEX
        static constexpr label nFaces() noexcept
        {
            return 6;
        }

        //- Return the model faces
        static const Foam::faceList& modelFaces();

        //- Return the model edges
        static const Foam::edgeList& modelEdges();


        //- Return i-th face
        inline Foam::face face(const label facei) const;

        //- Return i-th edge
        inline Foam::edge edge(const label edgei) const;

        //- Return i-th edge reversed
        inline Foam::edge reverseEdge(const label edgei) const;

        //- Return list of cell faces [6]
        Foam::faceList faces() const;

        //- Return list of cell edges [12]
        Foam::edgeList edges() const;

        //- Cell centre - uses simple average of points
        inline point centre(const UList<point>& meshPoints) const;

        //- The points corresponding to this shape
        inline pointField points(const UList<point>& meshPoints) const;


    // Operations

        //- Return HEX shape cell
        cellShape shape(const bool doCollapse = false) const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for hexCell
template<> struct is_contiguous<hexCell> : std::true_type {};

//- Contiguous label data for hexCell
template<> struct is_contiguous_label<hexCell> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hexCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
