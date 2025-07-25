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

#include "preservePatchesConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(preservePatches);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        preservePatches,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::preservePatches::preservePatches
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName),
    patches_(coeffDict_.get<wordRes>("patches"))
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep owner and (coupled) neighbour"
            << " of faces in patches " << patches_
            << " on same processor. This only makes sense for cyclics"
            << " and cyclicAMI." << endl;
    }
}


Foam::decompositionConstraints::preservePatches::preservePatches
(
    const UList<wordRe>& patches
)
:
    decompositionConstraint(dictionary(), typeName),
    patches_(patches)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep owner and (coupled) neighbour"
            << " of faces in patches " << patches_
            << " on same processor. This only makes sense for cyclics"
            << " and cyclicAMI." << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::preservePatches::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    blockedFace.resize(mesh.nFaces(), true);

    const labelList patchIDs(pbm.patchSet(patches_).sortedToc());

    label nUnblocked = 0;

    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, i)
        {
            const label meshFacei = pp.start() + i;

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


void Foam::decompositionConstraints::preservePatches::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // If the decomposition has not enforced the constraint, do it over here.

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const labelList patchIDs(pbm.patchSet(patches_).sortedToc());

    // Synchronise decomposition on patchIDs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nChanged = 0;

    do
    {
        // Extract min coupled boundary data
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList destProc;
        getMinBoundaryValue(mesh, decomposition, destProc);


        // Override (patchIDs only) if differing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nChanged = 0;
        for (const label patchi : patchIDs)
        {
            const polyPatch& pp = pbm[patchi];

            const labelUList& faceCells = pp.faceCells();

            forAll(faceCells, i)
            {
                const label bFacei = pp.offset()+i;
                if (destProc[bFacei] < decomposition[faceCells[i]])
                {
                    decomposition[faceCells[i]] = destProc[bFacei];
                    ++nChanged;
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
