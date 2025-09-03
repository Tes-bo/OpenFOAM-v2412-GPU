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

#include "kOmegaSSTDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow(Cd1_*this->r(this->nuEff(), magGradU), Cd2_));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDDES<BasicTurbulenceModel>::S2
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
        const volScalarField Ssigma(this->Ssigma(gradU));

        S2 -=
            fd(mag(gradU))
           *pos(this->lengthScaleRAS() - this->lengthScaleLES(CDES))
           *(S2 - sqr(Ssigma));
    }

    return tS2;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    const volScalarField lRAS(this->lengthScaleRAS());
    const volScalarField lLES(this->lengthScaleLES(CDES));
    const dimensionedScalar l0(dimLength, Zero);

    return max
    (
        lRAS - fd(magGradU)*max(lRAS - lLES, l0),
        dimensionedScalar("small", dimLength, SMALL)
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDDES<BasicTurbulenceModel>::GbyNu0
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
tmp<volScalarField::Internal> kOmegaSSTDDES<BasicTurbulenceModel>::GbyNu
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
kOmegaSSTDDES<BasicTurbulenceModel>::kOmegaSSTDDES
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
    kOmegaSSTDES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cd1_
    (
        this->useSigma_
          ? dimensioned<scalar>::getOrAddToDict
            (
                "Cd1Sigma",
                this->coeffDict_,
                22
            )
          : dimensioned<scalar>::getOrAddToDict
            (
                "Cd1",
                this->coeffDict_,
                20
            )
    ),
    Cd2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd2",
            this->coeffDict_,
            3
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTDDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTDES<BasicTurbulenceModel>::read())
    {
        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDDES<BasicTurbulenceModel>::fd() const
{
    return fd(mag(fvc::grad(this->U_)));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
