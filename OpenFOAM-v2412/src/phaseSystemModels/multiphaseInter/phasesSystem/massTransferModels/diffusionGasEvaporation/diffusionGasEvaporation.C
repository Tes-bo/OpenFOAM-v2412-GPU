/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "diffusionGasEvaporation.H"
#include "constants.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "fvcGrad.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::meltingEvaporationModels::
diffusionGasEvaporation<Thermo, OtherThermo>::updateInterface
(
    const volScalarField& T
)
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
        const label status = cutCell.calcSubCell(celli, isoAlpha_);
        interfaceArea_[celli] = 0;
        if (status == 0) // cell is cut
        {
            interfaceArea_[celli] = mag(cutCell.faceArea())/mesh.V()[celli];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::diffusionGasEvaporation<Thermo, OtherThermo>
::diffusionGasEvaporation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    saturationModelPtr_
    (
        saturationModel::New
        (
            dict.subDict("saturationPressure"),
            this->mesh_
        )
    ),
    isoAlpha_(dict.getOrDefault<scalar>("isoAlpha", 0.5)),
    C_("C", dimless, dict),
    Tactivate_("Tactivate", dimTemperature, 0, dict),
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
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::diffusionGasEvaporation<Thermo, OtherThermo>
::Kexp
(
    const volScalarField& T
)
{
    const fvMesh& mesh = this->mesh_;

    const word speciesName(IOobject::member(this->transferSpecie()));

    const typename OtherThermo::thermoType& vapourThermo =
        this->getLocalThermo
        (
            speciesName,
            this->toThermo_
        );

    const volScalarField& from = this->pair().from();
    const volScalarField& to = this->pair().to();

    const volScalarField& Yv = this->toThermo_.composition().Y(speciesName);

    updateInterface(T);

    auto tRhog = volScalarField::New
    (
        "tRhog",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimDensity, Zero)
    );
    auto& rhog = tRhog.ref();

    rhog = this->pair().to().rho();

    auto tDvg = volScalarField::New
    (
        "tDvg",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(sqr(dimLength)/dimTime, Zero)
    );
    auto& Dvg = tDvg.ref();

    Dvg = this->Dto(speciesName);

    tmp<volScalarField> tpSat = saturationModelPtr_->pSat(T);

    const volScalarField XvSat(tpSat()/this->toThermo_.p());

    const dimensionedScalar Wv("Wv", dimMass/dimMoles, vapourThermo.W());
    const volScalarField YvSat
    (
        XvSat
       *(
           Wv/(XvSat*Wv + (1-XvSat)*this->toThermo_.W())
        )
    );

    const volScalarField Ygm(max(from*YvSat + to*Yv, Zero));

    const multiphaseInterSystem& fluid = this->fluid();

    tmp<volVectorField> tnHatInt(fluid.nVolHatfv(to, from));

    const volScalarField gradYgm(fvc::grad(Ygm) & tnHatInt());

    mDotc_ =
       -pos(T - Tactivate_)
       *C_*rhog*Dvg*gradYgm*interfaceArea_
       /(1 - YvSat);

    if (debug && mesh.time().writeTime())
    {
        volScalarField pSat("pSat", saturationModelPtr_->pSat(T));
        pSat.write();

        volScalarField YvSat1("YvSat", YvSat);
        YvSat1.write();

        volScalarField YgmDebug("Ygm", Ygm);
        YgmDebug.write();

        volScalarField gradYgmD("gradYgm", gradYgm);
        gradYgmD.write();

        volVectorField nHatIntD("nHatInt", tnHatInt());
        nHatIntD.write();
    }

    return tmp<volScalarField>::New(mDotc_);
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::diffusionGasEvaporation<Thermo, OtherThermo>::
KSp
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::diffusionGasEvaporation<Thermo, OtherThermo>::
KSu
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


//************************************************************************ //
