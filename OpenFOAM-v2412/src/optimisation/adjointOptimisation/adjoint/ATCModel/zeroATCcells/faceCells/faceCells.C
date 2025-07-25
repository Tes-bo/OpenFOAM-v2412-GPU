/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "faceCells.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faceCells, 0);
addToRunTimeSelectionTable(zeroATCcells, faceCells, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faceCells::faceCells
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    zeroATCcells(mesh, dict)
{
    DynamicList<label> cellIDs;
    for (const fvPatch& patch : mesh_.boundary())
    {
        for (const word& patchType : zeroATCPatches_)
        {
            if (patch.type() == patchType)
            {
                cellIDs.push_back(patch.faceCells());
            }
        }
    }

    for (const label zoneID : zeroATCZones_)
    {
        if (zoneID != -1)
        {
            cellIDs.push_back(mesh_.cellZones()[zoneID]);
        }
    }
    zeroATCZones_.transfer(cellIDs);

    Info<< "Setting limiter on "
        << returnReduce(zeroATCcells_.size(), sumOp<label>())
        << " cells" << nl << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
