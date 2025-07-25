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

#include "hexCell.H"
#include "cellShape.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Warning:
// Ordering of faces needs to be the same for
// a hexahedral cell shape model and a hexCell

const Foam::label Foam::hexCell::modelFaces_[6][4] =
{
    {0, 4, 7, 3},  // x-min
    {1, 2, 6, 5},  // x-max
    {0, 1, 5, 4},  // y-min
    {3, 7, 6, 2},  // y-max
    {0, 3, 2, 1},  // z-min
    {4, 5, 6, 7}   // z-max
};


// Warning:
// Ordering of edges needs to be the same for
// a hexahedral cell shape model and a hexCell

const Foam::label Foam::hexCell::modelEdges_[12][2] =
{
    {0, 1},  // x-direction
    {3, 2},
    {7, 6},
    {4, 5},
    {0, 3},  // y-direction
    {1, 2},
    {5, 6},
    {4, 7},
    {0, 4},  // z-direction
    {1, 5},
    {2, 6},
    {3, 7}
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::faceList& Foam::hexCell::modelFaces()
{
    static std::unique_ptr<Foam::faceList> ptr(nullptr);

    if (!ptr)
    {
        ptr.reset(new Foam::faceList(hexCell::nFaces(), Foam::face(4)));

        label facei = 0;
        for (auto& f : *ptr)
        {
            f[0] = modelFaces_[facei][0];
            f[1] = modelFaces_[facei][1];
            f[2] = modelFaces_[facei][2];
            f[3] = modelFaces_[facei][3];
            ++facei;
        }
    }

    return *ptr;
}


const Foam::edgeList& Foam::hexCell::modelEdges()
{
    static std::unique_ptr<Foam::edgeList> ptr(nullptr);

    if (!ptr)
    {
        ptr.reset(new Foam::edgeList(hexCell::nEdges()));

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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::faceList Foam::hexCell::faces() const
{
    Foam::faceList theFaces(hexCell::nFaces(), Foam::face(4));

    label facei = 0;
    for (auto& f : theFaces)
    {
        copyFace(f, facei);
        ++facei;
    }

    return theFaces;
}


Foam::edgeList Foam::hexCell::edges() const
{
    Foam::edgeList theEdges(hexCell::nEdges());

    label edgei = 0;
    for (auto& e : theEdges)
    {
        e[0] = (*this)[modelEdges_[edgei][0]];
        e[1] = (*this)[modelEdges_[edgei][1]];
        ++edgei;
    }

    return theEdges;
}


Foam::cellShape Foam::hexCell::shape(const bool doCollapse) const
{
    static const cellModel* modelPtr(nullptr);

    if (!modelPtr)
    {
        modelPtr = cellModel::ptr(cellModel::HEX);
    }

    return cellShape(*modelPtr, *this, doCollapse);
}


// ************************************************************************* //
