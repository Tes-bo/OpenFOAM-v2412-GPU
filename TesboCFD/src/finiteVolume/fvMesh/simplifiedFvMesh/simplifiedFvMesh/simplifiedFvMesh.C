/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "simplifiedFvMesh.H"
#include "fvPatchField.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(simplifiedFvMesh, 0);
defineRunTimeSelectionTable(simplifiedFvMesh, time);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::simplifiedFvMesh::fvPatchFieldExists(const word& patchType)
{
    return
    (
        fvPatchField<scalar>::dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<vector>::dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<sphericalTensor>::
            dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<symmTensor>::
            dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<tensor>::dictionaryConstructorTablePtr_->found(patchType)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplifiedFvMesh::simplifiedFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour)
    )
{}


Foam::autoPtr<Foam::simplifiedFvMesh> Foam::simplifiedFvMesh::New
(
    const word& modelType,
    const Time& runTime
)
{
    Info<< "Selecting simplified mesh model " << modelType << endl;

    auto* ctorPtr = timeConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "simplified fvMesh",
            modelType,
            *timeConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<simplifiedFvMesh>(ctorPtr(runTime));
}


// ************************************************************************* //
