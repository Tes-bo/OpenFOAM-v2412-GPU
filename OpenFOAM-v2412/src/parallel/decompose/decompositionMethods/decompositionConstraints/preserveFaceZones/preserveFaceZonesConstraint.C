/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "preserveFaceZonesConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(preserveFaceZones);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        preserveFaceZones,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::preserveFaceZones::preserveFaceZones
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName),
    zones_(coeffDict_.get<wordRes>("zones"))
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : adding constraints to keep owner and neighbour"
            << " of faces in zones " << zones_
            << " on same processor" << endl;
    }
}


Foam::decompositionConstraints::preserveFaceZones::preserveFaceZones
(
    const UList<wordRe>& zones
)
:
    decompositionConstraint(dictionary(), typeName),
    zones_(zones)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : adding constraints to keep owner and neighbour"
            << " of faces in zones " << zones_
            << " on same processor" << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::preserveFaceZones::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    blockedFace.setSize(mesh.nFaces(), true);

    const faceZoneMesh& fZones = mesh.faceZones();

    const labelList zoneIDs(zones_.matching(fZones.names()));

    label nUnblocked = 0;

    for (const label zonei : zoneIDs)
    {
        const faceZone& fz = fZones[zonei];

        for (const label meshFacei : fz)
        {
            if (blockedFace[meshFacei])
            {
                blockedFace[meshFacei] = false;
                ++nUnblocked;
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        Info<< type() << " : unblocked "
            << returnReduce(nUnblocked, sumOp<label>()) << " faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::decompositionConstraints::preserveFaceZones::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // If the decomposition has not enforced the constraint do it over
    // here.

    label nChanged = 0;

    do
    {
        // Extract min coupled boundary data
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList destProc;
        getMinBoundaryValue(mesh, decomposition, destProc);


        // Override if differing
        // ~~~~~~~~~~~~~~~~~~~~~

        const faceZoneMesh& fZones = mesh.faceZones();

        const labelList zoneIDs(zones_.matching(fZones.names()));

        nChanged = 0;
        for (const label zonei : zoneIDs)
        {
            const faceZone& fz = fZones[zonei];

            for (const label facei : fz)
            {
                const label own = mesh.faceOwner()[facei];

                if (mesh.isInternalFace(facei))
                {
                    const label nei = mesh.faceNeighbour()[facei];
                    if (decomposition[nei] < decomposition[own])
                    {
                        decomposition[own] = decomposition[nei];
                        ++nChanged;
                    }
                }
                else
                {
                    const label bFaceI = facei-mesh.nInternalFaces();
                    if (destProc[bFaceI] < decomposition[own])
                    {
                        decomposition[own] = destProc[bFaceI];
                        ++nChanged;
                    }
                }
            }
        }

        if (decompositionConstraint::debug & 2)
        {
            Info<< type() << " : changed decomposition on "
                << returnReduce(nChanged, sumOp<label>()) << " cells" << endl;
        }

    } while (returnReduceOr(nChanged));
}


// ************************************************************************* //
