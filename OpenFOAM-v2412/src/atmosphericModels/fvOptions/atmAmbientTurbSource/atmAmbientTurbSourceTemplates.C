/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2023 OpenCFD Ltd
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

#include "atmAmbientTurbSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::atmAmbientTurbSource::atmAmbientTurbSourceEpsilon
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const tmp<volScalarField> tepsilon(turbPtr->epsilon());
    const volScalarField& epsilon = tepsilon();

    // (Heuristically derived from RS:Eq. 4, rhs-term:5)
    eqn +=
        fvm::Sp(alpha()*rho()*C2_*sqr(epsilonAmb_)/(kAmb_*epsilon()), epsilon);
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::atmAmbientTurbSource::atmAmbientTurbSourceOmega
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const tmp<volScalarField> tomega(turbPtr->omega());
    const volScalarField& omega = tomega();

    const volScalarField::Internal& beta =
        mesh_.lookupObjectRef<volScalarField::Internal>
        (
            IOobject::scopedName(turbPtr->type(), "beta")
        );

    // (RS:Eq. 4, rhs-term:5)
    eqn += fvm::Sp(alpha()*rho()*Cmu_*beta*sqr(omegaAmb_)/omega(), omega);
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::atmAmbientTurbSource::atmAmbientTurbSourceK
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const tmp<volScalarField> tk(turbPtr->k());
    const volScalarField& k = tk();

    if (isEpsilon_)
    {
        // (Heuristically derived from RS:Eq. 3, rhs-term:4)
        eqn += fvm::Sp(alpha()*rho()*epsilonAmb_/k(), k);
    }
    else
    {
        // (RS:Eq. 3, rhs-term:4)
        eqn += fvm::Sp(alpha()*rho()*Cmu_*omegaAmb_*kAmb_/k(), k);
    }
}


// ************************************************************************* //
