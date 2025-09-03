/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "SpalartAllmarasDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
const Foam::Enum
<
    typename Foam::LESModels::
        SpalartAllmarasDDES<BasicTurbulenceModel>::shieldingMode
>
Foam::LESModels::SpalartAllmarasDDES<BasicTurbulenceModel>::shieldingModeNames
({
    { shieldingMode::standard, "standard" },
    { shieldingMode::ZDES2020, "ZDES2020" },
});


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    const volScalarField r(this->r(this->nuEff(), magGradU, this->y_));

    tmp<volScalarField> tfd = 1 - tanh(pow(Cd1_*r, Cd2_));

    switch (shielding_)
    {
        case shieldingMode::standard:
        {
            return tfd;
        }
        case shieldingMode::ZDES2020:
        {
            auto maxEps = [](const volScalarField& fld, const scalar eps){
                return max(fld, dimensionedScalar(fld.dimensions(), eps));
            };

            volScalarField& fdStd = tfd.ref();
            const auto& nuTilda = this->nuTilda_;
            const volVectorField& n = wallDist::New(this->mesh_).n();

            const volScalarField GnuTilda
            (
                Cd3_*maxEps(fvc::grad(nuTilda) & n, Zero)
              / (maxEps(magGradU, SMALL)*this->kappa_*this->y_)
            );

            volScalarField fdGnuTilda(1 - tanh(pow(Cd1_*GnuTilda, Cd2_)));
            const volScalarField GOmega
            (
              - (fvc::grad(mag(fvc::curl(this->U_))) & n)
              * sqrt(nuTilda/maxEps(pow3(magGradU), SMALL))
            );
            const volScalarField alpha((7.0/6.0*Cd4_ - GOmega)/(Cd4_/6.0));
            const volScalarField fRGOmega
            (
                pos(Cd4_ - GOmega)
              + 1.0
               /(1 + exp(min(-6*alpha/max(1 - sqr(alpha), SMALL), scalar(50))))
               *pos(4*Cd4_/3.0 - GOmega)*pos(GOmega - Cd4_)
            );

            // Use more conservative fP2-function in case switch is true;
            // otherwise use simplified formulation
            if (usefP2_)
            {
                fdGnuTilda *=
                    (1.0 - tanh(pow(Cd1_*betaZDES_*r, Cd2_)))
		          / maxEps(fdStd, SMALL);
            }

            fdStd *= 1 - (1 - fdGnuTilda)*fRGOmega;
        }
    }

    return tfd;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU,
    const volScalarField& dTilda
) const
{
    if (this->useSigma_)
    {
        const volScalarField& lRAS(this->y_);
        const volScalarField fv2(this->fv2(chi, fv1));
        const volScalarField lLES(this->lengthScaleLES(chi, fv1));
        const volScalarField Omega(this->Omega(gradU));
        const volScalarField Ssigma(this->Ssigma(gradU));
        const volScalarField SsigmaDES
        (
            Omega - fd(mag(gradU))*pos(lRAS - lLES)*(Omega - Ssigma)
        );

        return
            max
            (
                SsigmaDES + fv2*this->nuTilda_/sqr(this->kappa_*dTilda),
                this->Cs_*SsigmaDES
            );
    }

    return
        SpalartAllmarasBase<DESModel<BasicTurbulenceModel>>::Stilda
        (
            chi,
            fv1,
            gradU,
            dTilda
        );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    const volScalarField& lRAS(this->y_);
    const volScalarField lLES(this->lengthScaleLES(chi, fv1));
    const dimensionedScalar l0(dimLength, Zero);

    return max
    (
        lRAS - fd(mag(gradU))*max(lRAS - lLES, l0),
        dimensionedScalar("small", dimLength, SMALL)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasDDES<BasicTurbulenceModel>::SpalartAllmarasDDES
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
    SpalartAllmarasDES<BasicTurbulenceModel>
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
    shielding_
    (
        shieldingModeNames.getOrDefault
        (
            "shielding",
            this->coeffDict_,
            shieldingMode::standard
        )
    ),
    Cd1_
    (
        this->useSigma_
          ? dimensioned<scalar>::getOrAddToDict
            (
                "Cd1Sigma",
                this->coeffDict_,
                10
            )
          : dimensioned<scalar>::getOrAddToDict
            (
                "Cd1",
                this->coeffDict_,
                8
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
    ),
    Cd3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd3",
            this->coeffDict_,
            25
        )
    ),
    Cd4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd4",
            this->coeffDict_,
            0.03
        )
    ),
    betaZDES_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaZDES",
            this->coeffDict_,
            2.5
        )
    ),
    usefP2_
    (
        Switch::getOrAddToDict
        (
            "usefP2",
            this->coeffDict_,
            false
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);

        switch (shielding_)
        {
            case shieldingMode::standard:
            {
                Info<< "shielding function: standard DDES "
                    <<  "(Spalart et al., 2006)"
                    << nl;
                break;
            }
            case shieldingMode::ZDES2020:
            {
                Info<< "shielding function: ZDES mode 2 (Deck & Renard, 2020)"
                    << nl;
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unrecognised 'shielding' option: "
                    << shieldingModeNames[shielding_]
                    << exit(FatalError);
            }
        }

        if (usefP2_)
        {
            Info<< "fP2 term: active" << nl;
        }
        else
        {
            Info<< "fP2 term: inactive" << nl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasDDES<BasicTurbulenceModel>::read()
{
    if (SpalartAllmarasDES<BasicTurbulenceModel>::read())
    {
        shieldingModeNames.readIfPresent
        (
            "shielding",
            this->coeffDict(),
            shielding_
        );

        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());
        Cd3_.readIfPresent(this->coeffDict());
        Cd4_.readIfPresent(this->coeffDict());
        betaZDES_.readIfPresent(this->coeffDict());
        usefP2_.readIfPresent("usefP2", this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::fd() const
{
    return fd(mag(fvc::grad(this->U_)));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
