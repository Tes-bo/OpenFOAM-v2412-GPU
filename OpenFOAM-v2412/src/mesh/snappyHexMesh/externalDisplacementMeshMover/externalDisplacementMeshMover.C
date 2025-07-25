/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "externalDisplacementMeshMover.H"
#include "mapPolyMesh.H"
#include "zeroFixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(externalDisplacementMeshMover, 0);
    defineRunTimeSelectionTable(externalDisplacementMeshMover, dictionary);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

Foam::labelList Foam::externalDisplacementMeshMover::getFixedValueBCs
(
    const pointVectorField& field
)
{
    DynamicList<label> adaptPatchIDs;

    forAll(field.boundaryField(), patchI)
    {
        const pointPatchField<vector>& patchField =
            field.boundaryField()[patchI];

        if (isA<valuePointPatchField<vector>>(patchField))
        {
            if (isA<zeroFixedValuePointPatchField<vector>>(patchField))
            {
                // Special condition of fixed boundary condition. Does not
                // get adapted
            }
            else
            {
                adaptPatchIDs.append(patchI);
            }
        }
    }

    return adaptPatchIDs;
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::externalDisplacementMeshMover::getPatch
(
    const polyMesh& mesh,
    const labelList& patchIDs
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Count faces.
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        nFaces += pp.size();
    }

    // Collect faces.
    labelList addressing(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        label meshFaceI = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFaceI++;
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>(mesh.faces(), addressing),
        mesh.points()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalDisplacementMeshMover::externalDisplacementMeshMover
(
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement,
    const bool dryRun
)
:
    baffles_(baffles),
    pointDisplacement_(pointDisplacement),
    dryRun_(dryRun)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::externalDisplacementMeshMover>
Foam::externalDisplacementMeshMover::New
(
    const word& type,
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement,
    const bool dryRun
)
{
    Info<< "Selecting externalDisplacementMeshMover " << type << endl;

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "externalDisplacementMeshMover",
            type,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<externalDisplacementMeshMover>
    (
        ctorPtr(dict, baffles, pointDisplacement, dryRun)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalDisplacementMeshMover::~externalDisplacementMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalDisplacementMeshMover::movePoints(const pointField&)
{
    // No local data to update
}


void Foam::externalDisplacementMeshMover::updateMesh(const mapPolyMesh& mpm)
{
    // Renumber baffles
    DynamicList<labelPair> newBaffles(baffles_.size());
    forAll(baffles_, i)
    {
        label f0 = mpm.reverseFaceMap()[baffles_[i].first()];
        label f1 = mpm.reverseFaceMap()[baffles_[i].second()];

        if (f0 >= 0 && f1 >= 0)
        {
            newBaffles.append(labelPair(f0, f1));
        }
    }
    newBaffles.shrink();
    baffles_.transfer(newBaffles);
}


// ************************************************************************* //
