/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::correctNut(S2);

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::r
(
    const volScalarField& nur,
    const volScalarField& magGradU
) const
{
    const dimensionedScalar eps(magGradU.dimensions(), SMALL);

    tmp<volScalarField> tr =
        min(nur/(max(magGradU, eps)*sqr(this->kappa_*this->y_)), scalar(10));

    tr.ref().boundaryFieldRef() == 0;

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::S2
(
    const volTensorField& gradU
) const
{
    tmp<volScalarField> tS2 =
        kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::S2(gradU);

    if (this->useSigma_)
    {
        volScalarField& S2 = tS2.ref();

        const volScalarField& k = this->k_;
        const volScalarField& omega = this->omega_;

        const volScalarField CDkOmega
        (
            (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
        );

        const volScalarField F1(this->F1(CDkOmega));

        const volScalarField CDES(this->CDES(F1));
        const volScalarField dTilda(this->dTilda(mag(gradU), CDES));
        const volScalarField lengthScaleRAS(this->lengthScaleRAS());
        const volScalarField Ssigma(this->Ssigma(gradU));

        S2 =
            pos(dTilda - lengthScaleRAS)*S2
          + (scalar(1) - pos(dTilda - lengthScaleRAS))*sqr(Ssigma);
    }

    return tS2;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    return min(lengthScaleLES(CDES), lengthScaleRAS());
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    volScalarField CDES(this->CDES(F1));
    return sqrt(this->k_())/dTilda(mag(gradU), CDES)()();
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicTurbulenceModel>::GbyNu0
(
    const volTensorField& gradU,
    const volScalarField& S2
) const
{
    if (this->useSigma_)
    {
        return S2();
    }

    return
        kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::GbyNu0(gradU, S2);
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicTurbulenceModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return GbyNu0; // Unlimited
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTDES<BasicTurbulenceModel>::kOmegaSSTDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTBase<DESModel<BasicTurbulenceModel>>
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

    useSigma_
    (
        Switch::getOrAddToDict
        (
            "useSigma",
            this->coeffDict_,
            false
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CDESkom_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CDESkom",
            this->coeffDict_,
            0.82
        )
    ),
    CDESkeps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CDESkeps",
            this->coeffDict_,
            0.60
        )
    )
{
    // Note: Ctrans coeff is model specific; for k-w = 60
    this->Ctrans_ =
        dimensioned<scalar>::getOrAddToDict("Ctrans", this->coeffDict_, 60.0);

    if (type == typeName)
    {
        WarningInFunction
            << "This model is not recommended and will be deprecated in future "
            << "releases. Please consider using the DDES/IDDES versions instead"
            << endl;

        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::read())
    {
        useSigma_.readIfPresent("useSigma", this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        CDESkom_.readIfPresent(this->coeffDict());
        CDESkeps_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
kOmegaSSTDES<BasicTurbulenceModel>::lengthScaleRAS() const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    return sqrt(k)/(this->betaStar_*omega);
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
kOmegaSSTDES<BasicTurbulenceModel>::lengthScaleLES
(
    const volScalarField& CDES
) const
{
    return CDES*this->delta();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volVectorField& U = this->U_;

    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    const volScalarField F1(this->F1(CDkOmega));

    return tmp<volScalarField>::New
    (
        IOobject::scopedName("DES", "LESRegion"),
        neg(dTilda(mag(fvc::grad(U)), CDES(F1)) - lengthScaleRAS())
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
