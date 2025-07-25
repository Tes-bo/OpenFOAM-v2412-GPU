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

#include "preserveBafflesConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "localPointRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(preserveBaffles);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        preserveBaffles,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::preserveBaffles::preserveBaffles
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : setting constraints to preserve baffles"
            //<< returnReduce(bafflePairs.size(), sumOp<label>())
            << endl;
    }
}


Foam::decompositionConstraints::preserveBaffles::preserveBaffles()
:
    decompositionConstraint(dictionary(), typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : setting constraints to preserve baffles"
            //<< returnReduce(bafflePairs.size(), sumOp<label>())
            << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::preserveBaffles::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    const labelPairList bafflePairs
    (
        localPointRegion::findDuplicateFacePairs(mesh)
    );

    if (decompositionConstraint::debug & 2)
    {
        Info<< type() << " : setting constraints to preserve "
            << returnReduce(bafflePairs.size(), sumOp<label>())
            << " baffles" << endl;
    }


    // Merge into explicitConnections
    {
        // Convert into face-to-face addressing
        labelList faceToFace(mesh.nFaces(), -1);
        for (const labelPair& p : explicitConnections)
        {
            faceToFace[p[0]] = p[1];
            faceToFace[p[1]] = p[0];
        }

        // Merge in bafflePairs
        for (const labelPair& p : bafflePairs)
        {
            if (faceToFace[p[0]] == -1 && faceToFace[p[1]] == -1)
            {
                faceToFace[p[0]] = p[1];
                faceToFace[p[1]] = p[0];
            }
            else if (labelPair::compare(p, labelPair(p[0], faceToFace[p[0]])))
            {
                // Connection already present
            }
            else
            {
                const label p0Slave = faceToFace[p[0]];
                const label p1Slave = faceToFace[p[1]];
                IOWarningInFunction(coeffDict_)
                    << "When adding baffle between faces "
                    << p[0] << " at " << mesh.faceCentres()[p[0]]
                    << " and "
                    << p[1] << " at " << mesh.faceCentres()[p[1]]
                    << " : face " << p[0] << " already is connected to face "
                    << p0Slave << " at " << mesh.faceCentres()[p0Slave]
                    << " and face " << p[1] << " already is connected to face "
                    << p1Slave << " at " << mesh.faceCentres()[p1Slave]
                    << endl;
            }
        }

        // Convert back into labelPairList
        label n = 0;
        forAll(faceToFace, faceI)
        {
            label otherFaceI = faceToFace[faceI];
            if (otherFaceI != -1 && faceI < otherFaceI)
            {
                // I am master of slave
                n++;
            }
        }
        explicitConnections.setSize(n);
        n = 0;
        forAll(faceToFace, faceI)
        {
            label otherFaceI = faceToFace[faceI];
            if (otherFaceI != -1 && faceI < otherFaceI)
            {
                explicitConnections[n++] = labelPair(faceI, otherFaceI);
            }
        }
    }

    // Make sure blockedFace is uptodate
    blockedFace.setSize(mesh.nFaces(), true);
    for (const labelPair& p : explicitConnections)
    {
        blockedFace[p.first()] = false;
        blockedFace[p.second()] = false;
    }
    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::decompositionConstraints::preserveBaffles::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    const labelPairList bafflePairs
    (
        localPointRegion::findDuplicateFacePairs(mesh)
    );

    label nChanged = 0;

    for (const labelPair& baffle : bafflePairs)
    {
        const label f0 = baffle.first();
        const label f1 = baffle.second();

        const label procI = decomposition[mesh.faceOwner()[f0]];

        if (mesh.isInternalFace(f0))
        {
            const label nei0 = mesh.faceNeighbour()[f0];
            if (decomposition[nei0] != procI)
            {
                decomposition[nei0] = procI;
                ++nChanged;
            }
        }

        const label own1 = mesh.faceOwner()[f1];
        if (decomposition[own1] != procI)
        {
            decomposition[own1] = procI;
            ++nChanged;
        }
        if (mesh.isInternalFace(f1))
        {
            const label nei1 = mesh.faceNeighbour()[f1];
            if (decomposition[nei1] != procI)
            {
                decomposition[nei1] = procI;
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        Info<< type() << " : changed decomposition on "
            << returnReduce(nChanged, sumOp<label>()) << " cells" << endl;
    }
}


// ************************************************************************* //
