/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "turbulentDispersionModel.H"
#include "phasePair.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(turbulentDispersionModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(turbulentDispersionModel, 0);
    defineRunTimeSelectionTable(turbulentDispersionModel, dictionary);
}

const Foam::dimensionSet Foam::turbulentDispersionModel::dimD(1, -1, -2, 0, 0);
const Foam::dimensionSet Foam::turbulentDispersionModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModel::turbulentDispersionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::turbulentDispersionModel>
Foam::turbulentDispersionModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting turbulentDispersionModel for "
        << pair << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "turbulentDispersionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return ctorPtr(dict, pair);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressibleTurbulenceModel&
Foam::turbulentDispersionModel::continuousTurbulence() const
{
    return
        pair_.phase1().mesh().lookupObject<phaseCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                pair_.continuous().name()
            )
        );
}


Foam::tmp<Foam::volVectorField>
Foam::turbulentDispersionModel::F() const
{
    return D()*fvc::grad(pair_.dispersed());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::turbulentDispersionModel::Ff() const
{
    return fvc::interpolate(D())*fvc::snGrad(pair_.dispersed());
}


// ************************************************************************* //
