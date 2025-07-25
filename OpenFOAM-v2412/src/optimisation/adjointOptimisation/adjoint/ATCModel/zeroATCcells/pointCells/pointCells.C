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

#include "pointCells.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointCells, 0);
addToRunTimeSelectionTable(zeroATCcells, pointCells, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointCells::pointCells
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    zeroATCcells(mesh, dict)
{
    boolList isZeroed(mesh_.nCells(), false);
    labelList zeroedIDs(mesh_.nCells(), -1);
    label i(0);
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        for (const word& patchType : zeroATCPatches_)
        {
            if (patch.type() == patchType)
            {
                const labelList& meshPoints =
                    mesh_.boundaryMesh()[patchI].meshPoints();

                for (const label pointI : meshPoints)
                {
                    const labelList& pointCells = mesh_.pointCells()[pointI];
                    for (const label cellI : pointCells)
                    {
                        if (!isZeroed[cellI])
                        {
                            zeroedIDs[i++] = cellI;
                            isZeroed[cellI] = true;
                        }
                    }
                }
            }
        }
    }
    forAll(zeroATCZones_, zI)
    {
        const label& zoneID = zeroATCZones_[zI];
        if (zoneID != -1)
        {
            const labelList& zoneCells = mesh_.cellZones()[zoneID];
            for (const label cellI : zoneCells)
            {
                if (!isZeroed[cellI])
                {
                    zeroedIDs[i++] = cellI;
                    isZeroed[cellI] = true;
                }
            }
        }
    }
    zeroedIDs.resize(i);
    zeroATCcells_ = std::move(zeroedIDs);

    Info<< "Zeroing ATC on "
        << returnReduce(zeroATCcells_.size(), sumOp<label>())
        << " cells" << nl << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
