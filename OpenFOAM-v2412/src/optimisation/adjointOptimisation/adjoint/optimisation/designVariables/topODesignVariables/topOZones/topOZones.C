/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "topOZones.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topOZones, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::topOZones::getZoneIDs
(
    const dictionary& dict,
    const word& zoneGroup
)
{
    wordList zoneNames(dict.getOrDefault<wordList>(zoneGroup, wordList(0)));
    labelList IDs(zoneNames.size(), -1);
    forAll(zoneNames, zI)
    {
        IDs[zI] = mesh_.cellZones().findZoneID(zoneNames[zI]);
    }
    return IDs;
}


void Foam::topOZones::addIOcellsZone()
{
    // Gather IO cells
    DynamicList<label> IOcells(0);
    for (const fvPatch& patch : mesh_.boundary())
    {
        if (patch.type() == "patch")
        {
            IOcells.append(patch.faceCells());
        }
    }

    // Add zone to cellZoneMesh and populate it
    cellZoneMesh& cellZones = const_cast<fvMesh&>(mesh_).cellZones();
    cellZone& IOcellsZone = cellZones("IOcells", true);
    IOcellsZone = IOcells;
    IOcellsID_ = cellZones.size() - 1;

    // Write as set
    cellSet IOcellList(mesh_, "IOcellList", IOcells);
    IOcellList.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topOZones::topOZones
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    fixedPorousIDs_(getZoneIDs(dict, "fixedPorousZones")),
    fixedPorousValues_
    (
        dict.getOrDefault<scalarList>
        (
            "fixedPorousValues", scalarList(fixedPorousIDs_.size(), 1.)
        )
    ),
    fixedZeroPorousIDs_(getZoneIDs(dict, "fixedZeroPorousZones")),
    adjointPorousIDs_(getZoneIDs(dict, "adjointPorousZones")),
    IOcellsID_(-1),
    betaMaxPtr_(betaMax::New(mesh, dict))
{
    addIOcellsZone();
    if (fixedPorousValues_.size() != fixedPorousIDs_.size())
    {
        FatalErrorInFunction
            << "Number of fixedPorousValues and fixedPorousZones don't agree!"
            << "\nSize  of fixedPorousIDs is " << fixedPorousIDs_.size()
            << " and size of fixedPorousValues is " << fixedPorousValues_.size()
            << endl << endl
            << exit(FatalError);
    }
}


// ************************************************************************* //
