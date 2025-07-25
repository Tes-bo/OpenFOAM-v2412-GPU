/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "cellPointWeightWallModified.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::cellPointWeightWallModified::onWall
(
    const polyMesh& mesh,
    const label facei
)
{
    if (facei >= 0)
    {
        const polyBoundaryMesh& bm = mesh.boundaryMesh();
        const label patchi = bm.whichPatch(facei);

        if (patchi != -1 && isA<wallPolyPatch>(bm[patchi]))
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointWeightWallModified::cellPointWeightWallModified
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label facei
)
:
    cellPointWeight(mesh, position, celli, facei)
{
    if (facei >= 0 && cellPointWeightWallModified::onWall(mesh, facei))
    {
        // Apply cell centre value for wall faces
        weights_[0] = 1;
        weights_[1] = 0;
        weights_[2] = 0;
        weights_[3] = 0;
    }
}


// ************************************************************************* //
