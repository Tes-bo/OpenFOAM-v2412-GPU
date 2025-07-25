/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "laminarModel.H"
#include "Stokes.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::laminarModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::laminarModel<BasicTurbulenceModel>::laminarModel
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

    laminarDict_(this->subOrEmptyDict("laminar")),
    printCoeffs_(laminarDict_.getOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(laminarDict_.optionalSubDict(type + "Coeffs"))
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::laminarModel<BasicTurbulenceModel>>
Foam::laminarModel<BasicTurbulenceModel>::New
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

    const dictionary* dictptr = modelDict.findDict("laminar");

    if (dictptr)
    {
        const dictionary& dict = *dictptr;

        const word modelType
        (
            // laminarModel -> model (after v2006)
            dict.getCompat<word>("model", {{"laminarModel", -2006}})
        );

        Info<< "Selecting laminar stress model " << modelType << endl;

        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "laminar model",
                modelType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<laminarModel>
        (
            ctorPtr
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport, propertiesName)
        );
    }
    else
    {
        Info<< "Selecting laminar stress model "
            << laminarModels::Stokes<BasicTurbulenceModel>::typeName << endl;

        return autoPtr<laminarModel>
        (
            new laminarModels::Stokes<BasicTurbulenceModel>
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport,
                propertiesName
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::laminarModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        laminarDict_ <<= this->subDict("laminar");

        coeffDict_ <<= laminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::nut() const
{
    return volScalarField::New
    (
        IOobject::groupName("nut", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(dimViscosity, Zero)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminarModel<BasicTurbulenceModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>::New(this->mesh_.boundary()[patchi].size(), Zero);
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::nuEff() const
{
    return volScalarField::New
    (
        IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->nu()
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminarModel<BasicTurbulenceModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::k() const
{
    return volScalarField::New
    (
        IOobject::groupName("k", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions())/dimTime, Zero)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::omega() const
{
    return volScalarField::New
    (
        IOobject::groupName("omega", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::laminarModel<BasicTurbulenceModel>::R() const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("R", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedSymmTensor(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicTurbulenceModel>
void Foam::laminarModel<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
