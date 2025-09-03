/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline Foam::pointField& Foam::edgeMesh::storedPoints() noexcept
{
    return points_;
}


inline Foam::edgeList& Foam::edgeMesh::storedEdges() noexcept
{
    return edges_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::edgeMesh::edgeMesh()
:
    fileFormats::edgeMeshFormatsCore(),
    points_(),
    edges_(),
    pointEdgesPtr_(nullptr)
{}


inline Foam::edgeMesh::edgeMesh(const edgeMesh& em)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(em.points_),
    edges_(em.edges_),
    pointEdgesPtr_(nullptr)
{}


inline Foam::edgeMesh::edgeMesh(edgeMesh&& em)
:
    edgeMesh()
{
    transfer(em);
}


inline Foam::edgeMesh::edgeMesh
(
    const pointField& points,
    const edgeList& edges
)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(points),
    edges_(edges),
    pointEdgesPtr_(nullptr)
{}


inline Foam::edgeMesh::edgeMesh
(
    pointField&& points,
    edgeList&& edges
)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(std::move(points)),
    edges_(std::move(edges)),
    pointEdgesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::pointField& Foam::edgeMesh::points() const noexcept
{
    return points_;
}


inline const Foam::edgeList& Foam::edgeMesh::edges() const noexcept
{
    return edges_;
}


inline const Foam::labelListList& Foam::edgeMesh::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        calcPointEdges();
    }
    return *pointEdgesPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::edgeMesh::operator=(const edgeMesh& rhs)
{
    points_ = rhs.points_;
    edges_ = rhs.edges_;
    pointEdgesPtr_.reset(nullptr);
}


inline void Foam::edgeMesh::operator=(edgeMesh&& rhs)
{
    transfer(rhs);
}


// ************************************************************************* //
