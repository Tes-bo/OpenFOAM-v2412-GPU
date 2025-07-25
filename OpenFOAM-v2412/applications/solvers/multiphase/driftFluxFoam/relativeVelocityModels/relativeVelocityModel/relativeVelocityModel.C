/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "relativeVelocityModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

Foam::wordList Foam::relativeVelocityModel::UdmPatchFieldTypes() const
{
    const volVectorField& U = mixture_.U();

    wordList UdmTypes
    (
        U.boundaryField().size(),
        fvPatchFieldBase::calculatedType()
    );

    forAll(U.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
        )
        {
            UdmTypes[i] = fixedValueFvPatchVectorField::typeName;
        }
    }

    return UdmTypes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    mixture_(mixture),
    alphac_(mixture.alpha2()),
    alphad_(mixture.alpha1()),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),

    Udm_
    (
        IOobject
        (
            "Udm",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector(dimVelocity, Zero),
        UdmPatchFieldTypes()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
{
    const word modelType(dict.get<word>(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "relative velocity",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return
        autoPtr<relativeVelocityModel>
        (
            ctorPtr
            (
                dict.optionalSubDict(modelType + "Coeffs"),
                mixture
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::relativeVelocityModel::rho() const
{
    return alphac_*rhoc_ + alphad_*rhod_;
}


Foam::tmp<Foam::volSymmTensorField> Foam::relativeVelocityModel::tauDm() const
{
    volScalarField betac(alphac_*rhoc_);
    volScalarField betad(alphad_*rhod_);

    // Calculate the relative velocity of the continuous phase w.r.t the mean
    volVectorField Ucm(betad*Udm_/betac);

    return volSymmTensorField::New
    (
        "tauDm",
        IOobject::NO_REGISTER,
        betad*sqr(Udm_) + betac*sqr(Ucm)
    );
}


// ************************************************************************* //
