/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "pointSmoother.H"
#include "pointConstraints.H"
#include "polyMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointSmoother, 0);
    defineRunTimeSelectionTable(pointSmoother, dictionary);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

bool Foam::pointSmoother::isInternalOrProcessorFace(const label faceI) const
{
    if (mesh().isInternalFace(faceI))
    {
        return true;
    }

    if
    (
        processorPatchIDs_
        [
            mesh().boundaryMesh().patchID()
            [
                faceI - mesh().nInternalFaces()
            ]
        ]
    )
    {
        return true;
    }

    return false;
}


Foam::autoPtr<Foam::PackedBoolList> Foam::pointSmoother::pointsToMove
(
    const labelList& facesToMove,
    const bool moveInternalFaces
) const
{
    autoPtr<PackedBoolList> markerPtr
    (
        new PackedBoolList(mesh().nPoints(), false)
    );

    PackedBoolList& marker(markerPtr());

    forAll(facesToMove, faceToMoveI)
    {
        const label faceI(facesToMove[faceToMoveI]);

        if (moveInternalFaces || !isInternalOrProcessorFace(faceI))
        {
            const face& fPoints(mesh().faces()[faceI]);

            forAll(fPoints, fPointI)
            {
                const label pointI(fPoints[fPointI]);

                marker[pointI] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        marker,
        orEqOp<unsigned int>(),
        0U
    );

    return markerPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointSmoother::pointSmoother
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh)
{
    for (const auto& pp : mesh.boundaryMesh())
    {
        if (isA<processorPolyPatch>(pp))
        {
            processorPatchIDs_.insert(pp.index());
        }
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pointSmoother>
Foam::pointSmoother::New
(
    const word& pointSmootherType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    Info<< "Selecting pointSmoother type " << pointSmootherType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(pointSmootherType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("pointSmoother::New")
            << "Unknown " << typeName << " type "
            << pointSmootherType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, dict);
}


Foam::autoPtr<Foam::pointSmoother>
Foam::pointSmoother::New
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    word pointSmootherType(dict.lookup(typeName));

    return New(pointSmootherType, mesh, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointSmoother::~pointSmoother()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointSmoother::update
(
    const labelList& facesToMove,
    const pointField& oldPoints,
    const pointField& currentPoints,
    const pointField& faceCentres,
    const vectorField& faceAreas,
    const pointField& cellCentres,
    const scalarField& cellVolumes,
    pointVectorField& pointDisplacement,
    const bool correctBCs
) const
{
    // Apply point smoothing (without boundary conditions!)
    calculate
    (
        facesToMove,
        oldPoints,
        currentPoints,
        faceCentres,
        faceAreas,
        cellCentres,
        cellVolumes,
        pointDisplacement.ref()
    );

    // Note: do not want to apply boundary conditions whilst smoothing so
    //       disable for now        
    //// Take over boundary values
    //for (auto& ppf : pointDisplacement.boundaryFieldRef())
    //{
    //    ppf.operator==(ppf.patchInternalField());
    //}

    if (correctBCs)
    {
        // Apply (multi-)patch constraints
        pointConstraints::New
        (
            pointDisplacement().mesh()
        ).constrainDisplacement(pointDisplacement);
    }
}


Foam::tmp<Foam::scalarField> Foam::pointSmoother::faceQuality
(
    const pointField& points,
    const pointField& faceCentres,
    const vectorField& faceAreas,
    const pointField& cellCentres,
    const scalarField& cellVolumes
) const
{
    // See e.g.
    // - src/functionObjects/field/stabilityBlendingFactor/
    //      stabilityBlendingFactor.H
    // - maybe add function to motionSmootherCheck to return value instead
    //   of yes/no for invalid?
    // - for now just look at non-ortho

    tmp<scalarField> tortho
    (
        polyMeshTools::faceOrthogonality
        (
            mesh(),
            faceAreas,
            cellCentres
        )
    );

    //return max(tortho, scalar(0.0));
    return tortho;
}


Foam::tmp<Foam::scalarField> Foam::pointSmoother::cellQuality
(
    const pointField& points,
    const pointField& faceCentres,
    const vectorField& faceAreas,
    const pointField& cellCentres,
    const scalarField& cellVolumes
) const
{
    // Get face-based non-ortho
    tmp<scalarField> tfaceOrtho
    (
        faceQuality
        (
            points,
            faceCentres,
            faceAreas,
            cellCentres,
            cellVolumes
        )
    );
    const auto& faceOrtho = tfaceOrtho();

    // Min over cells
    tmp<scalarField> tortho(new scalarField(mesh().nCells(), GREAT));
    auto& ortho = tortho.ref();

    const auto& own = mesh().faceOwner();
    const auto& nei = mesh().faceNeighbour();

    forAll(own, facei)
    {
        auto& o = ortho[own[facei]];
        o = min(o, faceOrtho[facei]);
    }
    for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
    {
        auto& o = ortho[nei[facei]];
        o = min(o, faceOrtho[facei]);
    }

    return tortho;
}


// ************************************************************************* //
