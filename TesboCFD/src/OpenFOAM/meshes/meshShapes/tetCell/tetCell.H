/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::tetCell

Description
    A tetrahedral cell primitive.

    It is important that the face/edge ordering is identical for
    a tetrahedron class,
    a tetrahedron cell shape model and a tetCell

SourceFiles
    tetCell.C
    tetCellI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_tetCell_H
#define Foam_tetCell_H

#include "FixedList.H"
#include "triFace.H"
#include "faceList.H"
#include "edgeList.H"
#include "pointField.H"
#include "tetrahedron.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellShape;

/*---------------------------------------------------------------------------*\
                           Class tetCell Declaration
\*---------------------------------------------------------------------------*/

class tetCell
:
    public FixedList<label, 4>
{
    // Static Data Members

        //- The model faces for TET
        static const label modelFaces_[4][3];

        //- The model edges for TET
        static const label modelEdges_[6][2];


public:

    // Generated Methods

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Default construct, with invalid point labels (-1)
        inline tetCell();

        //- Construct from four point labels
        inline tetCell
        (
            const label p0,
            const label p1,
            const label p2,
            const label p3
        );

        //- Construct from an initializer list of four point labels
        inline explicit tetCell(std::initializer_list<label> list);

        //- Construct from FixedList of four point labels
        inline tetCell(const FixedList<label, 4>& list);

        //- Copy construct from a subset of point labels
        inline tetCell
        (
            const labelUList& list,
            const FixedList<label, 4>& indices
        );

        //- Copy construct from a subset of point labels
        template<unsigned AnyNum>
        inline tetCell
        (
            const FixedList<label, AnyNum>& list,
            const FixedList<label, 4>& indices
        );

        //- Construct from Istream
        inline explicit tetCell(Istream& is);


    // Member Functions

        //- Number of points for TET
        static constexpr label nPoints() noexcept
        {
            return 4;
        }

        //- Number of edges for TET
        static constexpr label nEdges() noexcept
        {
            return 6;
        }

        //- Number of faces for TET
        static constexpr label nFaces() noexcept
        {
            return 4;
        }

        //- Return the model faces
        static const Foam::faceList& modelFaces();

        //- Return the model edges
        static const Foam::edgeList& modelEdges();


        //- Return i-th face
        inline Foam::triFace face(const label facei) const;

        //- Return first face adjacent to the given edge
        inline label edgeFace(const label edgei) const;

        //- Return face adjacent to the given face sharing the same edge
        inline label edgeAdjacentFace
        (
            const label edgei,
            const label facei
        ) const;

        //- Return i-th edge
        inline Foam::edge edge(const label edgei) const;

        //- Return i-th edge reversed
        inline Foam::edge reverseEdge(const label edgei) const;

        //- The points corresponding to this shape
        inline pointField points(const UList<point>& meshPoints) const;


    // Operations

        //- Return TET shape cell
        cellShape shape() const;

        //- Return the tetrahedron
        inline tetPointRef tet(const UList<point>& meshPoints) const;


    // Housekeeping

        //- Identical to edge()
        Foam::edge tetEdge(label edgei) const { return this->edge(edgei); }

        //- Identical to shape()
        cellShape tetCellShape() const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for tetCell
template<> struct is_contiguous<tetCell> : std::true_type {};

//- Contiguous label data for tetCell
template<> struct is_contiguous_label<tetCell> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
