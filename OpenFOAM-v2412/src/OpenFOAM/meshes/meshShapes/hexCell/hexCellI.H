/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::hexCell::copyFace(Foam::face& f, const label facei) const
{
    f[0] = (*this)[modelFaces_[facei][0]];
    f[1] = (*this)[modelFaces_[facei][1]];
    f[2] = (*this)[modelFaces_[facei][2]];
    f[3] = (*this)[modelFaces_[facei][3]];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::hexCell::hexCell()
:
    FixedList<label, 8>(-1)
{}


inline Foam::hexCell::hexCell(std::initializer_list<label> list)
:
    FixedList<label, 8>(list)
{}


inline Foam::hexCell::hexCell(const FixedList<label, 8>& list)
:
    FixedList<label, 8>(list)
{}


inline Foam::hexCell::hexCell
(
    const labelUList& list,
    const FixedList<label, 8>& indices
)
:
    FixedList<label, 8>(list, indices)
{}


template<unsigned AnyNum>
inline Foam::hexCell::hexCell
(
    const FixedList<label, AnyNum>& list,
    const FixedList<label, 8>& indices
)
:
    FixedList<label, 8>(list, indices)
{}


inline Foam::hexCell::hexCell(Istream& is)
:
    FixedList<label, 8>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::face Foam::hexCell::face(const label facei) const
{
    #ifdef FULLDEBUG
    if (facei < 0 || facei >= hexCell::nFaces())
    {
        FatalErrorInFunction
            << "Face index (" << facei << ") out of range 0..5\n"
            << abort(FatalError);
    }
    #endif

    Foam::face f(4);
    copyFace(f, facei);

    return f;
}


inline Foam::edge Foam::hexCell::edge(const label edgei) const
{
    #ifdef FULLDEBUG
    if (edgei < 0 || edgei >= hexCell::nEdges())
    {
        FatalErrorInFunction
            << "Edge index (" << edgei << ") out of range 0..11\n"
            << abort(FatalError);
    }
    #endif

    return Foam::edge
    (
        (*this)[modelEdges_[edgei][0]],
        (*this)[modelEdges_[edgei][1]]
    );
}


inline Foam::edge Foam::hexCell::reverseEdge(const label edgei) const
{
    // Reverse edge. Using a copy is cheaper than inplace flip
    return this->edge(edgei).reverseEdge();
}


inline Foam::point Foam::hexCell::centre
(
    const UList<point>& meshPoints
) const
{
    // Simple estimate of cell centre by averaging cell points
    point cEst = Zero;
    int npts = 0;
    for (const label pointi : *this)
    {
        if (pointi >= 0)
        {
            cEst += meshPoints[pointi];
            ++npts;
        }
    }

    return (npts > 1 ? (cEst/scalar(npts)) : cEst);
}


inline Foam::pointField Foam::hexCell::points
(
    const UList<point>& meshPoints
) const
{
    return pointField(List<point>(meshPoints, *this));
}


// ************************************************************************* //
