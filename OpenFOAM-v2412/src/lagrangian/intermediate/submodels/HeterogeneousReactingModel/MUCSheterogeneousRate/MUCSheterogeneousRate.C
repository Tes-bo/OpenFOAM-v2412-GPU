/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "MUCSheterogeneousRate.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MUCSheterogeneousRate<CloudType>::MUCSheterogeneousRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    HeterogeneousReactingModel<CloudType>(dict, owner, typeName),
    D12_(this->coeffDict().getScalar("D12")),
    epsilon_(this->coeffDict().getScalar("epsilon")),
    gamma_(this->coeffDict().getScalar("gamma")),
    sigma_(this->coeffDict().getScalar("sigma")),
    E_(this->coeffDict().getScalar("E")),
    A_(this->coeffDict().getScalar("A")),
    Aeff_(this->coeffDict().getScalar("Aeff")),
    Ea_(this->coeffDict().getScalar("Ea")),
    nuFuel_(this->coeffDict().getScalar("nuFuel")),
    nuOx_(this->coeffDict().getScalar("nuOx")),
    nuProd_(this->coeffDict().getScalar("nuProd")),
    O2GlobalId_(owner.composition().carrierId("O2")),
    FuelLocalId_(-1),
    ProdLocalId_(-1),
    WO2_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    FuelLocalId_ =
        owner.composition().localId
        (
            idSolid,
            this->coeffDict().getWord("fuel")
        );

    ProdLocalId_ =
        owner.composition().localId
        (
            idSolid,
            this->coeffDict().getWord("product")
        );

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().W(O2GlobalId_);
}


template<class CloudType>
Foam::MUCSheterogeneousRate<CloudType>::MUCSheterogeneousRate
(
    const MUCSheterogeneousRate<CloudType>& srm
)
:
    HeterogeneousReactingModel<CloudType>(srm),
    D12_(srm.D12_),
    epsilon_(srm.epsilon_),
    gamma_(srm.gamma_),
    sigma_(srm.sigma_),
    E_(srm.E_),
    A_(srm.A_),
    Aeff_(srm.Aeff_),
    Ea_(srm.Ea_),
    nuFuel_(srm.nuFuel_),
    nuOx_(srm.nuOx_),
    nuProd_(srm.nuProd_),
    O2GlobalId_(srm.O2GlobalId_),
    FuelLocalId_(srm.FuelLocalId_),
    ProdLocalId_(srm.ProdLocalId_),
    WO2_(srm.WO2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::MUCSheterogeneousRate<CloudType>::calculate
(
    const scalar dt,
    const scalar Re,
    const scalar nu,
    const label celli,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    const scalarField& YSolid,
    scalarField& F,
    const scalar N,
    scalar& NCpW,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier
) const
{
    // Fraction of remaining combustible material
    const scalar fComb = YSolid[FuelLocalId_];

    // Surface combustion until combustible fraction is consumed
    if (fComb < SMALL)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();
    const auto& composition = this->owner().composition();

    const scalar WFuel = composition.solids().properties()[FuelLocalId_].W();
    const scalar WProd = composition.solids().properties()[ProdLocalId_].W();

    // O2 concentration [Kmol/m3]
    const scalar Cb =
        thermo.carrier().Y(O2GlobalId_)[celli]*rhoc/WO2_;

    if (Cb < SMALL)
    {
        return 0.0;
    }

    // Reaction constant
    const scalar k = A_*exp(-Ea_/(RR*T));

    // Effective diffussivity
    const scalar Deff = D12_*epsilon_/gamma_;

     // Schmidt number
    const scalar Sc = nu/(D12_ + ROOTVSMALL);

    // Mass transfer coefficient [m/s]
    const scalar alpha =
        (2.0 + 0.6*Foam::sqrt(Re)*cbrt(Sc))*D12_/(d + ROOTVSMALL);

    const scalar r = d/2;

    const scalar f = F[FuelLocalId_];

    const scalar rhof = composition.solids().properties()[FuelLocalId_].rho();

    const scalar deltaRho0 = (nuOx_/nuFuel_)*rhof/WFuel;

    // Progress variable rate
    const scalar dfdt =
        Aeff_*(Cb/deltaRho0)
       /(
           r/3/alpha
         + sqr(r)*(1/cbrt(1-f)-1)/3/Deff
         - (1/sqr(cbrt(1-f)))*r/k/sigma_/E_/3
        );

    // Update new progress variable
    F[FuelLocalId_] += dfdt*dt;

    // Interface radius
    const scalar ri = r*cbrt(1-f);

    // Interface radius rate
    //const scalar dridt = -dfdt*(r/3)*pow(1-f, -2/3);
    const scalar dridt = -dfdt*(pow3(r)/3)/sqr(ri);

    // O2 flux [Kmol/sec]
    const scalar q02 = deltaRho0*4*constant::mathematical::pi*sqr(ri)*dridt;

    // Calculate the number of molar units reacted [Kmol]
    const scalar dOmega = q02*dt;

    // Heat of Reaction
    const scalar Hc =
        composition.solids().properties()[ProdLocalId_].Hf()
      - composition.solids().properties()[FuelLocalId_].Hf();

    //Stoichiometric mass ratio for fuel
    const scalar sFuel = nuFuel_/(nuOx_);

    //Stoichiometric mass ratio for product
    const scalar sProd = nuProd_/(nuOx_);

    // Add to carrier phase mass transfer [Kg]
    dMassSRCarrier[O2GlobalId_] += dOmega*WO2_;

    // Remove to particle mass transfer
    dMassSolid[FuelLocalId_] -= dOmega*WFuel*sFuel;

    // Add to particle product
    dMassSolid[ProdLocalId_] += dOmega*WProd*sProd;

    if (debug)
    {
        Pout<< "mass    = " << mass << nl
            << "fComb   = " << fComb << nl
            << "dfdt    = " << dfdt << nl
            << "F       = " << F[FuelLocalId_] << nl
            << "ri      = " << ri << nl
            << "dridt   = " << dridt << nl
            << "q02     = " << q02 << nl
            << "dOmega  = " << dOmega << nl
            << "Hr      = " << dOmega*WFuel*sFuel*Hc << endl;
    }

    // Heat of reaction [J]
    return -dOmega*WFuel*sFuel*Hc;
}


template<class CloudType>
Foam::label Foam::MUCSheterogeneousRate<CloudType>::nReactions() const
{
    return 1;
}


// ************************************************************************* //
