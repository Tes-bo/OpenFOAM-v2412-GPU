/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
#include "polyMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cuttingSurfaceBase::checkOverlap
(
    const word callerName,
    const boundBox& meshBounds,
    const boundBox& userBounds
)
{
    // User bounding-box does not overlap with (global) mesh!
    if (userBounds.good() && !userBounds.overlaps(meshBounds))
    {
        WarningInFunction
            << nl << callerName
            << " : Bounds " << userBounds
            << " do not overlap the mesh bounding box " << meshBounds
            << nl << endl;
    }
}


Foam::bitSet Foam::cuttingSurfaceBase::cellSelection
(
    const polyMesh& mesh,
    const boundBox& userBounds,
    const wordRes& zoneNames,
    boundBox& meshBounds
)
{
    bitSet cellsToSelect;

    // Zones requested and in use?
    const bool hasZones =
        returnReduceAnd(-1 != mesh.cellZones().findIndex(zoneNames));

    if (hasZones)
    {
        cellsToSelect = mesh.cellZones().selection(zoneNames);
    }


    // Subset the zoned cells with the userBounds.
    // For a full mesh, use the bounds to define the cell selection.

    // If there are zones cells, use them to build the effective mesh
    // bound box.
    // Note that for convenience we use cell centres here instead of
    // cell points, since it will only be used for checking.


    meshBounds = mesh.bounds();  // Use the regular mesh bounding box

    const auto& cellCentres = static_cast<const fvMesh&>(mesh).C();

    if (userBounds.empty())
    {
        // No bounds restriction, but may need effective mesh
        // bounding-box for later checks

        if (hasZones)
        {
            meshBounds.clear();

            for (const label celli : cellsToSelect)
            {
                const point& cc = cellCentres[celli];

                meshBounds.add(cc);
            }

            meshBounds.reduce();
        }
    }
    else if (hasZones)
    {
        // Subset zoned cells with the user bounding-box

        for (const label celli : cellsToSelect)
        {
            const point& cc = cellCentres[celli];

            meshBounds.add(cc);

            if (!userBounds.contains(cc))
            {
                cellsToSelect.unset(celli);
            }
        }

        meshBounds.reduce();
    }
    else
    {
        // Create cell selection from user bounding-box

        const label len = mesh.nCells();

        cellsToSelect.resize(len);

        for (label celli=0; celli < len; ++celli)
        {
            const point& cc = cellCentres[celli];

            if (userBounds.contains(cc))
            {
                cellsToSelect.set(celli);
            }
        }
    }

    return cellsToSelect;
}


Foam::bitSet Foam::cuttingSurfaceBase::cellSelection
(
    const polyMesh& mesh,
    const boundBox& userBounds,
    const wordRes& zoneNames,
    const word callerName,
    const bool warn
)
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
