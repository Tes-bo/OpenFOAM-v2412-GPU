/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "viewFactorModel.H"
#include "viewFactorHottel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::VF::viewFactorModel> Foam::VF::viewFactorModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // Intercept 2D cases
    if (mesh.nSolutionD() == 2)
    {
        Info<< "Selecting " << typeName << ": " << viewFactorHottel::typeName
            << " for 2D cases" << nl << endl;
        return autoPtr<viewFactorModel>(new viewFactorHottel(mesh, dict));
    }


    // 3D cases - use run-time selection

    const word modelType(dict.get<word>("viewFactorModel"));

    Info<< "Selecting " << typeName << ": " << modelType << endl;

    auto* ctorPtr = meshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            typeName,
            modelType,
            *meshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<viewFactorModel>(ctorPtr(mesh, dict));
}


// ************************************************************************* //
