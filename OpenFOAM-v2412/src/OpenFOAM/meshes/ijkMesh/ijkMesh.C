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

#include "ijkMesh.H"
#include "hexCell.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::hexCell
Foam::ijkMesh::vertLabels(const label i, const label j, const label k) const
{
    hexCell verts;

    verts[0] = pointLabel(i,   j,   k);
    verts[1] = pointLabel(i+1, j,   k);
    verts[2] = pointLabel(i+1, j+1, k);
    verts[3] = pointLabel(i,   j+1, k);
    verts[4] = pointLabel(i,   j,   k+1);
    verts[5] = pointLabel(i+1, j,   k+1);
    verts[6] = pointLabel(i+1, j+1, k+1);
    verts[7] = pointLabel(i,   j+1, k+1);

    return verts;
}


Foam::hexCell Foam::ijkMesh::vertLabels(const labelVector& ijk) const
{
    return vertLabels(ijk.x(), ijk.y(), ijk.z());
}


// ************************************************************************* //
