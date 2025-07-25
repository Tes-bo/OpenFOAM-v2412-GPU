/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
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

#include "DESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
DESModel<BasicTurbulenceModel>::DESModel
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
    LESeddyViscosity<BasicTurbulenceModel>
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
    Ctrans_(dimless, Zero)
{
    // Note: Ctrans is model-specific and initialised in derived classes
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool DESModel<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ctrans_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> DESModel<BasicTurbulenceModel>::Ssigma
(
    const volTensorField& gradU,
    const dimensionedScalar& coeff
)
{
    // Limiter
    const dimensionedScalar eps0(dimless, SMALL);
    const dimensionedScalar eps2(dimless/sqr(dimTime), SMALL);
    const dimensionedScalar eps4(dimless/pow4(dimTime), SMALL);
    const dimensionedScalar max2(dimless/sqr(dimTime), GREAT);
    const dimensionedTensor maxTen2
    (
        dimless/sqr(dimTime),
        tensor::max
    );
    const dimensionedTensor minTen2
    (
        dimless/sqr(dimTime),
        tensor::min
    );

    const volTensorField G(max(min(gradU.T() & gradU, maxTen2), minTen2));

    // Tensor invariants
    const volScalarField I1(tr(G));
    const volScalarField I2(0.5*(sqr(I1) - tr(G & G)));
    tmp<volScalarField> tI3 = det(G);

    const volScalarField alpha1(max(sqr(I1)/9.0 - I2/3.0, eps4));

    tmp<volScalarField> talpha2 =
        pow3(min(I1, max2))/27.0 - I1*I2/6.0 + 0.5*tI3;

    const volScalarField alpha3
    (
        1.0/3.0
       *acos
        (
            max
            (
                scalar(-1) + eps0,
                min(scalar(1) - eps0, talpha2/pow(alpha1, 3.0/2.0))
            )
        )
    );

    const scalar piBy3 = constant::mathematical::pi/3.0;
    const volScalarField sigma1
    (
        sqrt(max(I1/3.0 + 2.0*sqrt(alpha1)*cos(alpha3), eps2))
    );
    const volScalarField sigma2
    (
        sqrt(max(I1/3.0 - 2.0*sqrt(alpha1)*cos(piBy3 + alpha3), eps2))
    );
    const volScalarField sigma3
    (
        sqrt(max(I1/3.0 - 2.0*sqrt(alpha1)*cos(piBy3 - alpha3), eps2))
    );

    return
        coeff
       *sigma3
       *(sigma1 - sigma2)
       *(sigma2 - sigma3)
       /max(sqr(sigma1), eps2);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> DESModel<BasicTurbulenceModel>::Ssigma
(
    const volTensorField& gradU
) const
{
    return Ssigma(gradU, Ctrans_);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> DESModel<BasicTurbulenceModel>::fd() const
{
    return volScalarField::New
    (
        "fd",
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(dimless, Zero)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
