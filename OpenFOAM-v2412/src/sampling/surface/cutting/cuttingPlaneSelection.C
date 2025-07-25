/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
#include "polyMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cuttingPlane::checkOverlap
(
    const word callerName,
    const boundBox& meshBounds,
    const boundBox& userBounds
) const
{
    cuttingSurfaceBase::checkOverlap(callerName, meshBounds, userBounds);

    const plane& pln = *this;

    // Plane does not intersect the user bounding box
    if (userBounds.good() && !userBounds.intersects(pln))
    {
        WarningInFunction
            << nl << callerName
            << " : Plane "<< pln << " does not intersect the bounds "
            << userBounds
            << nl << endl;
    }

    // Plane does not intersect the (global) mesh!
    if (!meshBounds.intersects(pln))
    {
        WarningInFunction
            << nl << callerName
            << " : Plane "<< pln << " does not intersect the mesh bounds "
            << meshBounds
            << nl << endl;
    }
}


Foam::bitSet Foam::cuttingPlane::cellSelection
(
    const polyMesh& mesh,
    const boundBox& userBounds,
    const wordRes& zoneNames,
    const word callerName,
    const bool warn
) const
{
    boundBox meshBounds;

    bitSet cellsToSelect =
        cuttingSurfaceBase::cellSelection
        (
            mesh, userBounds, zoneNames, meshBounds
        );

    if (warn)
    {
        checkOverlap(callerName, meshBounds, userBounds);
    }

    return cellsToSelect;
}


// ************************************************************************* //
