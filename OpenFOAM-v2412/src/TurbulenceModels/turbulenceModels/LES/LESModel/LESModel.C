/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "LESModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::LESModel<BasicTurbulenceModel>::LESModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    LESDict_(this->subOrEmptyDict("LES")),
    turbulence_(LESDict_.getOrDefault<Switch>("turbulence", true)),
    printCoeffs_(LESDict_.getOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(LESDict_.optionalSubDict(type + "Coeffs")),
    Ce_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ce",
            LESDict_,
            1.048
        )
    ),
    kMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kMin",
            LESDict_,
            sqr(dimVelocity),
            SMALL
        )
    ),
    epsilonMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "epsilonMin",
            LESDict_,
            kMin_.dimensions()/dimTime,
            SMALL
        )
    ),
    omegaMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaMin",
            LESDict_,
            dimless/dimTime,
            SMALL
        )
    ),
    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            LESDict_
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::LESModel<BasicTurbulenceModel>>
Foam::LESModel<BasicTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    const IOdictionary modelDict
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const dictionary& dict = modelDict.subDict("LES");

    const word modelType
    (
        // "LESModel" for v2006 and earlier
        dict.getCompat<word>("model", {{"LESModel", -2006}})
    );

    Info<< "Selecting LES turbulence model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "LES model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<LESModel>
    (
        ctorPtr(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::LESModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        LESDict_ <<= this->subDict("LES");
        LESDict_.readEntry("turbulence", turbulence_);

        coeffDict_ <<= LESDict_.optionalSubDict(type() + "Coeffs");

        delta_().read(LESDict_);

        Ce_.readIfPresent(LESDict_);

        kMin_.readIfPresent(LESDict_);

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::LESModel<BasicTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        this->Ce()*pow(this->k(), 1.5)/this->delta()
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::LESModel<BasicTurbulenceModel>::omega() const
{
    const scalar betaStar = 0.09;
    const dimensionedScalar k0(sqr(dimLength/dimTime), SMALL);

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        this->epsilon()/(betaStar*(this->k() + k0))
    );
}


template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::correct()
{
    delta_().correct();
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
