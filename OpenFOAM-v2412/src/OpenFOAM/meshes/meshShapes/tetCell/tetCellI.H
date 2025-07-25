/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::tetCell::tetCell()
:
    FixedList<label, 4>(-1)
{}


inline Foam::tetCell::tetCell
(
    const label p0,
    const label p1,
    const label p2,
    const label p3
)
{
    operator[](0) = p0;
    operator[](1) = p1;
    operator[](2) = p2;
    operator[](3) = p3;
}


inline Foam::tetCell::tetCell(std::initializer_list<label> list)
:
    FixedList<label, 4>(list)
{}


inline Foam::tetCell::tetCell(const FixedList<label, 4>& list)
:
    FixedList<label, 4>(list)
{}


inline Foam::tetCell::tetCell
(
    const labelUList& list,
    const FixedList<label, 4>& indices
)
:
    FixedList<label, 4>(list, indices)
{}


template<unsigned AnyNum>
inline Foam::tetCell::tetCell
(
    const FixedList<label, AnyNum>& list,
    const FixedList<label, 4>& indices
)
:
    FixedList<label, 4>(list, indices)
{}


inline Foam::tetCell::tetCell(Istream& is)
:
    FixedList<label, 4>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::triFace Foam::tetCell::face(const label facei) const
{
    #ifdef FULLDEBUG
    if (facei < 0 || facei >= tetCell::nFaces())
    {
        FatalErrorInFunction
            << "Face index (" << facei << ") out of range 0..3\n"
            << abort(FatalError);
    }
    #endif

    return Foam::triFace
    (
        (*this)[modelFaces_[facei][0]],
        (*this)[modelFaces_[facei][1]],
        (*this)[modelFaces_[facei][2]]
    );
}


inline Foam::label Foam::tetCell::edgeFace(const label edgei) const
{
    // Warning. Ordering of faces needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell
    static const label edgeFaces[6] = {2, 3, 1, 0, 0, 1};

    #ifdef FULLDEBUG
    if (edgei < 0 || edgei >= tetCell::nEdges())
    {
        FatalErrorInFunction
            << "Edge index (" << edgei << ") out of range 0..5\n"
            << abort(FatalError);
    }
    #endif

    return edgeFaces[edgei];
}


inline Foam::label Foam::tetCell::edgeAdjacentFace
(
    const label edgei,
    const label facei
) const
{
    // Warning. Ordering of faces needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell
    static const label adjacentFace[6][4] =
    {
        {-1, -1, 3, 2},
        {-1, 3, -1, 1},
        {-1, 2, 1, -1},
        {2, -1, 0, -1},
        {3, -1, -1, 0},
        {1, 0, -1, -1}
    };

    #ifdef FULLDEBUG
    if (facei < 0 || facei >= tetCell::nFaces())
    {
        FatalErrorInFunction
            << "Face index (" << facei << ") out of range 0..3\n"
            << abort(FatalError);
    }

    if (edgei < 0 || edgei >= tetCell::nEdges())
    {
        FatalErrorInFunction
            << "Edge index (" << edgei << ") out of range 0..5\n"
            << abort(FatalError);
    }
    #endif

    return adjacentFace[edgei][facei];
}


inline Foam::edge Foam::tetCell::edge(const label edgei) const
{
    #ifdef FULLDEBUG
    if (edgei < 0 || edgei >= tetCell::nEdges())
    {
        FatalErrorInFunction
            << "Edge index (" << edgei << ") out of range 0..5\n"
            << abort(FatalError);
    }
    #endif

    return Foam::edge
    (
        (*this)[modelEdges_[edgei][0]],
        (*this)[modelEdges_[edgei][1]]
    );
}


inline Foam::edge Foam::tetCell::reverseEdge(const label edgei) const
{
    // Reverse edge. Using a copy is cheaper than inplace flip
    return this->edge(edgei).reverseEdge();
}


inline Foam::pointField Foam::tetCell::points
(
    const UList<point>& meshPoints
) const
{
    return pointField(List<point>(meshPoints, *this));
}


inline Foam::tetPointRef Foam::tetCell::tet
(
    const UList<point>& meshPoints
) const
{
    return tetPointRef(meshPoints, *this);
}


// ************************************************************************* //
