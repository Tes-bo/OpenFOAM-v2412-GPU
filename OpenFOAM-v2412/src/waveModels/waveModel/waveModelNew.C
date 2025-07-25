/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "waveModel.H"
#include "fvMesh.H"

Foam::autoPtr<Foam::waveModel> Foam::waveModel::New
(
    const word& dictName,
    const fvMesh& mesh,
    const polyPatch& patch
)
{
    IOdictionary waveDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    word modelType("none");
    dictionary patchDict;
    if (waveDict.found(patch.name()))
    {
        patchDict = waveDict.subDict(patch.name());
        modelType = patchDict.get<word>("waveModel");
    }
    else
    {
        FatalIOErrorInFunction(waveDict)
            << "Dictionary entry for patch " << patch.name() << " not found"
            << exit(FatalIOError);
    }

    Info<< "Selecting waveModel " << modelType << endl;

    auto* ctorPtr = patchConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            waveDict,
            "waveModel",
            modelType,
            *patchConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<waveModel>(ctorPtr(patchDict, mesh, patch));
}


Foam::tmp<Foam::waveModel> Foam::waveModel::lookupOrCreate
(
    const polyPatch& patch,
    const fvMesh& mesh,
    const word& waveDictName
)
{
    const word modelName = waveModel::modelName(patch.name());

    waveModel* modelPtr = mesh.getObjectPtr<waveModel>(modelName);

    if (!modelPtr)
    {
        modelPtr = waveModel::New(waveDictName, mesh, patch).ptr();
        modelPtr->store();
        modelPtr->info(Info);
    }

    return *modelPtr;
}


// ************************************************************************* //
