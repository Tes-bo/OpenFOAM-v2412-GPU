/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "surfaceTensionModel.H"
#include "constantSurfaceTension.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceTensionModel> Foam::surfaceTensionModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const dictionary* sigmaDictPtr = dict.findDict("sigma");

    if (sigmaDictPtr)
    {
        const dictionary& sigmaDict = *sigmaDictPtr;

        const word modelType(sigmaDict.get<word>("type"));

        Info<< "Selecting surfaceTensionModel " << modelType << endl;

        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "surfaceTensionModel",
                modelType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return ctorPtr(sigmaDict, mesh);
    }

    return autoPtr<surfaceTensionModel>
    (
        new surfaceTensionModels::constant(dict, mesh)
    );
}


// ************************************************************************* //
