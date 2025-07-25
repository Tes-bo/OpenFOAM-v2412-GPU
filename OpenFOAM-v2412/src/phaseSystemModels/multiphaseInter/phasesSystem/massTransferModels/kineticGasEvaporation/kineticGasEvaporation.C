/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "kineticGasEvaporation.H"
#include "constants.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "fvcSmooth.H"

using namespace Foam::constant;


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::updateInterface(const volScalarField& T)
{
    const fvMesh& mesh = this->mesh_;

    const volScalarField& alpha = this->pair().from();

    scalarField ap
    (
        volPointInterpolation::New(mesh).interpolate(alpha)
    );

    cutCellIso cutCell(mesh, ap);

    forAll(interfaceArea_, celli)
    {
        label status = cutCell.calcSubCell(celli, isoAlpha_);
        interfaceArea_[celli] = 0;
        if (status == 0) // cell is cut
        {
            interfaceArea_[celli] =
                mag(cutCell.faceArea())/mesh.V()[celli];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::kineticGasEvaporation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C", dimless, dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    Mv_("Mv", dimMass/dimMoles, -1, dict),
    interfaceArea_
    (
        IOobject
        (
            "interfaceArea",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    htc_
    (
        IOobject
        (
            "htc",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimMass/dimArea/dimTemperature/dimTime, Zero)
    ),
    mDotc_
    (
        IOobject
        (
            "mDotc",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    isoAlpha_(dict.getOrDefault<scalar>("isoAlpha", 0.5))
{
    word speciesName = IOobject::member(this->transferSpecie());

    // Get the "to" thermo
    const typename OtherThermo::thermoType& toThermo =
        this->getLocalThermo
        (
            speciesName,
            this->toThermo_
        );

    // Convert from g/mol to Kg/mol
    Mv_.value() = toThermo.W()*1e-3;

    if (Mv_.value() == -1)
    {
        FatalErrorInFunction
            << " Please provide the molar weight (Mv) of vapour [g/mol] "
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Kexp(const volScalarField& T)
{

    const fvMesh& mesh = this->mesh_;

    const dimensionedScalar HerztKnudsConst
    (
        sqrt
        (
            2.0*mathematical::pi
          * pow3(Tactivate_)
          * constant::physicoChemical::R/Mv_
        )
    );

    word speciesName = IOobject::member(this->transferSpecie());
    tmp<volScalarField> L = mag(this->L(speciesName, T));

    updateInterface(T);

    auto tRhov = volScalarField::New
    (
        "tRhov",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimDensity, Zero)
    );
    auto& rhov = tRhov.ref();

    auto tdeltaT = volScalarField::New
    (
        "tdeltaT",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimTemperature, Zero)
    );
    auto& deltaT = tdeltaT.ref();

    dimensionedScalar T0("T0", dimTemperature, Zero);

    if (sign(C_.value()) > 0)
    {
        rhov = this->pair().to().rho();
        deltaT = max(T - Tactivate_, T0);
    }
    else
    {
        rhov = this->pair().from().rho();
        deltaT = max(Tactivate_ - T, T0);
    }

    htc_ = 2*mag(C_)/(2-mag(C_))*(L()*rhov/HerztKnudsConst);

    mDotc_ = htc_*deltaT*interfaceArea_;

    return tmp<volScalarField>::New(mDotc_);
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>::KSp
(
    label variable,
    const volScalarField& refValue
)
{
    if (this->modelVariable_ == variable)
    {
        const volScalarField coeff(htc_*interfaceArea_);

        if (sign(C_.value()) > 0)
        {
            return(coeff*pos(refValue - Tactivate_));
        }
        else
        {
            return(coeff*pos(Tactivate_ - refValue));
        }
    }

    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>::KSu
(
    label variable,
    const volScalarField& refValue
)
{
    if (this->modelVariable_ == variable)
    {
        const volScalarField coeff(htc_*interfaceArea_*Tactivate_);

        if (sign(C_.value()) > 0)
        {
            return(-coeff*pos(refValue - Tactivate_));
        }
        else
        {
            return(coeff*pos(Tactivate_ - refValue));
        }
    }

    return nullptr;
}


// ************************************************************************* //
