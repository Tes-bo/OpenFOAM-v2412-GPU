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

#include "cuttingSurfaceBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cuttingSurfaceBase::debug
(
    Foam::debug::debugSwitch("cuttingSurfaceBase", 0)
);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cuttingSurfaceBase::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const bitSet& cellIdLabels
)
{
    bitSet subsetCells(cellIdLabels);

    performCut(mesh, triangulate, std::move(subsetCells));
}


void Foam::cuttingSurfaceBase::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const labelUList& cellIdLabels
)
{
    bitSet subsetCells;

    if (notNull(cellIdLabels))
    {
        // Pre-populate with restriction
        subsetCells.resize(mesh.nCells());
        subsetCells.set(cellIdLabels);
    }

    performCut(mesh, triangulate, std::move(subsetCells));
}


void Foam::cuttingSurfaceBase::remapFaces(const labelUList& faceMap)
{
    if (!faceMap.empty())
    {
        Mesh::remapFaces(faceMap);

        List<label> remappedCells(faceMap.size());
        forAll(faceMap, facei)
        {
            remappedCells[facei] = meshCells_[faceMap[facei]];
        }
        meshCells_.transfer(remappedCells);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cuttingSurfaceBase::operator=(const cuttingSurfaceBase& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    static_cast<Mesh&>(*this) = rhs;
    meshCells_ = rhs.meshCells();
}


// ************************************************************************* //
