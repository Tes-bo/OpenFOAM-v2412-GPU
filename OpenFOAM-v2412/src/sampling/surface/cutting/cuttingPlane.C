/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "cuttingPlane.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttingPlane::cuttingPlane(const plane& pln)
:
    plane(pln)
{}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    const bitSet& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    bitSet&& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    const labelUList& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cuttingPlane::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    bitSet&& cellIdLabels
)
{
    const plane& pln = *this;
    const pointField& pts = mesh.points();

    Mesh::clear();
    meshCells_.clear();

    // Pre-populate with restriction
    bitSet cellCuts(std::move(cellIdLabels));

    if (cellCuts.size())
    {
        cellCuts.resize(mesh.nCells());
    }

    // For each mesh point, the encoded side (0,1,2) of the plane
    PackedList<2> sides;

    // Determine cells that are (likely) cut
    // - some ambiguity when plane is exactly between cells
    const label nFaceCuts = calcCellCuts(mesh, sides, cellCuts);


    // Walk cell cuts to create faces

    // Action #1:
    // - Orient edge so it points in the positive normal direction.
    // - Edge/plane intersection when the sign changes
    const auto edgeOrientIntersect =
        [=](edge& e) -> bool
        {
            if (sides[e.last()] < sides[e.first()])
            {
                e.flip();
            }

            return sides[e.first()] != sides[e.last()];
        };

    // Action #2:
    // - The edge intersection alpha
    const auto edgeAlphaIntersect =
        [=](const edge& e) -> scalar
        {
            return pln.lineIntersect(e.line(pts));
        };

    walkCellCuts
    (
        mesh,
        cellCuts,
        edgeOrientIntersect,
        edgeAlphaIntersect,
        triangulate,
        nFaceCuts
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cuttingPlane::operator=(const cuttingPlane& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    static_cast<Mesh&>(*this) = rhs;
    static_cast<plane&>(*this) = rhs;
    meshCells_ = rhs.meshCells();
}


// ************************************************************************* //
