/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "dynamicTopODesignVariables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicTopODesignVariables, 0);
    addToRunTimeSelectionTable
    (
        designVariables,
        dynamicTopODesignVariables,
        designVariables
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::dynamicTopODesignVariables::setActiveDesignVariables
(
    const label fluidID,
    const bool activeIO
)
{
    cellZoneMesh& cellZones = mesh_.cellZones();
    marchCells_.addFixedCells(cellZones, zones_.fixedPorousZoneIDs());
    marchCells_.addFixedCells(cellZones, zones_.fixedZeroPorousZoneIDs());
    if (!activeIO)
    {
        marchCells_.addFixedCells(zones_.IOCells());
    }

    // Evolve the active design variables as many times as the passed
    // optimisation cycles - 1 (the last evolution corresponds to the design
    // variables of the next optimisation cycle).
    // Ensure at least one evolution, to kick-start the optimisation
    evolvedCount_ = max(evolvedCount_ - 1, 1);
    marchCells_.update(evolvedCount_);
    activeDesignVariables_ = marchCells_.getActiveCells();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicTopODesignVariables::dynamicTopODesignVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    dynamicTopODesignVariables(mesh, dict, mesh.nCells())
{}


Foam::dynamicTopODesignVariables::dynamicTopODesignVariables
(
    fvMesh& mesh,
    const dictionary& dict,
    const label size
)
:
    topODesignVariables(mesh, dict, size),
    marchCells_(mesh, dict.subDict("marchingCoeffs")),
    evolvedCount_(localIOdictionary::getOrDefault<label>("evolvedCount", 0))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicTopODesignVariables>
Foam::dynamicTopODesignVariables::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    return autoPtr<dynamicTopODesignVariables>::New(mesh, dict);
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::dynamicTopODesignVariables::evolveNumber()
{
    marchCells_.update();
    activeDesignVariables_ = marchCells_.getActiveCells();
    ++evolvedCount_;
}


bool Foam::dynamicTopODesignVariables::writeData(Ostream& os) const
{
    os.writeEntry("evolvedCount", evolvedCount_);
    return topODesignVariables::writeData(os);
}


// ************************************************************************* //
