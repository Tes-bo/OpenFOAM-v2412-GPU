/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "zoneMotion.H"
#include "syncTools.H"
#include "bitSet.H"
#include "cellSet.H"
#include "cellZoneMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneMotion::zoneMotion
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    pointIDs_(),
    moveAllCells_(true)
{
    // Specified cellSet?
    word cellSetName;

    if
    (
        dict.readIfPresent("cellSet", cellSetName)
     && cellSetName == "none"  // Compat: ignore 'none' placeholder
    )
    {
        cellSetName.clear();
    }

    labelList cellIDs;
    if (!cellSetName.empty())
    {
        Info<< "Applying motion to cellSet: " << cellSetName << endl;

        cellIDs = cellSet(mesh, cellSetName).toc();
    }


    // Specified cellZone(s) ?
    wordRe cellZoneName;

    if
    (
        dict.readIfPresent("cellZone", cellZoneName)
     && cellZoneName == "none"  // Compat: ignore 'none' placeholder
    )
    {
        cellZoneName.clear();
    }

    labelList zoneIDs;
    if (!cellZoneName.empty())
    {
        Info<< "Applying motion to cellZone: " << cellZoneName << endl;

        // Also handles groups, multiple zones (as wordRe match) ...
        zoneIDs = mesh.cellZones().indices(cellZoneName);

        if (zoneIDs.empty())
        {
            FatalIOErrorInFunction(dict)
                << "No matching cellZones: " << cellZoneName << nl
                << "    Valid zones : "
                << flatOutput(mesh.cellZones().names()) << nl
                << "    Valid groups: "
                << flatOutput(mesh.cellZones().groupNames())
                << nl
                << exit(FatalIOError);
        }
    }

    if (!cellSetName.empty() || !cellZoneName.empty())
    {
        bitSet movePts(mesh.nPoints());

        // Markup points associated with cell zone(s)
        for (const label zoneID : zoneIDs)
        {
            for (const label celli : mesh.cellZones()[zoneID])
            {
                for (const label facei : mesh.cells()[celli])
                {
                    movePts.set(mesh.faces()[facei]);
                }
            }
        }

        // Markup points associated with cellSet
        for (const label celli : cellIDs)
        {
            for (const label facei : mesh.cells()[celli])
            {
                movePts.set(mesh.faces()[facei]);
            }
        }

        syncTools::syncPointList(mesh, movePts, orEqOp<unsigned int>(), 0u);

        pointIDs_ = movePts.sortedToc();
    }


    // No cell points selected (as set or zones) => move all points

    moveAllCells_ = returnReduceAnd(pointIDs_.empty());

    if (moveAllCells_)
    {
        Info<< "Applying motion to entire mesh" << endl;
    }
}


// ************************************************************************* //
