/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "ErgunWenYuDragForce.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ErgunWenYuDragForce<CloudType>::CdRe(const scalar alphacRe)
const
{
    // (ZZB:Eq. 14, GLSLR:Table 3)
    if (alphacRe < 1000.0)
    {
        return 24.0*(1.0 + 0.15*pow(alphacRe, 0.687));
    }

    return 0.44*alphacRe;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ErgunWenYuDragForce<CloudType>::ErgunWenYuDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    alphac_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            this->coeffs().getWord("alphac")
        )
    )
{}


template<class CloudType>
Foam::ErgunWenYuDragForce<CloudType>::ErgunWenYuDragForce
(
    const ErgunWenYuDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    alphac_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            this->coeffs().getWord("alphac")
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::ErgunWenYuDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar alphac = alphac_[p.cell()];

    if (alphac < 0.8)
    {
        return forceSuSp
        (
            Zero,
            (mass/p.rho())
           *(150.0*(1.0 - alphac)/alphac + 1.75*Re)*muc/(alphac*sqr(p.d()))
        );
    }
    else
    {
        return forceSuSp
        (
            Zero,
            (mass/p.rho())
           *0.75*CdRe(alphac*Re)*muc*pow(alphac, -2.65)/(alphac*sqr(p.d()))
        );
    }
}


// ************************************************************************* //
