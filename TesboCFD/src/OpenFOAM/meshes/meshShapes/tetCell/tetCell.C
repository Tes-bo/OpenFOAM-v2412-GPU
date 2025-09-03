/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "tetCell.H"
#include "cellShape.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Warning.
// Ordering of faces needs to be the same for
// a tetrahedron class, a tetrahedron cell shape model and a tetCell

const Foam::label Foam::tetCell::modelFaces_[4][3] =
{
    {1, 2, 3},
    {0, 3, 2},
    {0, 1, 3},
    {0, 2, 1},
};


// Warning.
// Ordering of edges needs to be the same for
// a tetrahedron class, a tetrahedron cell shape model and a tetCell

const Foam::label Foam::tetCell::modelEdges_[6][2] =
{
    {0, 1},
    {0, 2},
    {0, 3},
    {3, 1},
    {1, 2},
    {3, 2}
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::faceList& Foam::tetCell::modelFaces()
{
    static std::unique_ptr<Foam::faceList> ptr(nullptr);

    if (!ptr)
    {
        ptr.reset(new Foam::faceList(tetCell::nFaces(), Foam::face(3)));

        label facei = 0;
        for (auto& f : *ptr)
        {
            f[0] = modelFaces_[facei][0];
            f[1] = modelFaces_[facei][1];
            f[2] = modelFaces_[facei][2];
            ++facei;
        }
    }

    return *ptr;
}


const Foam::edgeList& Foam::tetCell::modelEdges()
{
    static std::unique_ptr<Foam::edgeList> ptr(nullptr);

    if (!ptr)
    {
        ptr.reset(new Foam::edgeList(tetCell::nEdges()));

        label edgei = 0;
        for (auto& e : *ptr)
        {
            e[0] = modelEdges_[edgei][0];
            e[1] = modelEdges_[edgei][1];
            ++edgei;
        }
    }

    return *ptr;
}


/// Foam::faceList Foam::tetCell::faces() const
/// {
///     Foam::faceList theFaces(tetCell::nFaces(), Foam::face(3));
///
///     label facei = 0;
///     for (auto& f : theFaces)
///     {
///         f[0] = (*this)[modelFaces_[facei][0]];
///         f[1] = (*this)[modelFaces_[facei][1]];
///         f[2] = (*this)[modelFaces_[facei][2]];
///         ++facei;
///     }
///
///     return theFaces;
/// }
///
///
/// Foam::edgeList Foam::tetCell::edges() const
/// {
///     Foam::edgeList theEdges(tetCell::nEdges());
///
///     label edgei = 0;
///     for (auto& e : theEdges)
///     {
///         e[0] = (*this)[modelEdges_[edgei][0]];
///         e[1] = (*this)[modelEdges_[edgei][1]];
///         ++edgei;
///     }
///
///     return theEdges;
/// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::cellShape Foam::tetCell::shape() const
{
    static const cellModel* modelPtr(nullptr);

    if (!modelPtr)
    {
        modelPtr = cellModel::ptr(cellModel::TET);
    }

    return cellShape(*modelPtr, *this);
}


Foam::cellShape Foam::tetCell::tetCellShape() const
{
    return this->shape();
}


// ************************************************************************* //
