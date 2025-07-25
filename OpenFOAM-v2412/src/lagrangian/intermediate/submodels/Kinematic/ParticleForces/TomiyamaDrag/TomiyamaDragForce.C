/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "TomiyamaDragForce.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::Enum<typename Foam::TomiyamaDragForce<CloudType>::contaminationType>
Foam::TomiyamaDragForce<CloudType>::contaminationTypeNames
{
    { contaminationType::PURE, "pure" },
    { contaminationType::SLIGHT, "slight" },
    { contaminationType::FULL, "full" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::TomiyamaDragForce<CloudType>::CdRe(const scalar Re) const
{
    const scalar f = 1 + 0.15*pow(Re, 0.687);

    switch (contaminationType_)
    {
        case contaminationType::PURE:
        {
            // Eq. 31 pure system
            return min(16*f, 48);
        }
        case contaminationType::SLIGHT:
        {
            // Eq. 32 slightly contaminated system
            return min(24*f, 72);
        }
        case contaminationType::FULL:
        {
            // Eq. 33 fully contaminated system
            return 24*f;
        }
        default:
        {}
    }

    return Zero;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TomiyamaDragForce<CloudType>::TomiyamaDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    sigma_(this->coeffs().getScalar("sigma")),
    contaminationType_
    (
        contaminationTypeNames.get("contamination", this->coeffs())
    )
{}


template<class CloudType>
Foam::TomiyamaDragForce<CloudType>::TomiyamaDragForce
(
    const TomiyamaDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    sigma_(df.sigma_),
    contaminationType_(df.contaminationType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::TomiyamaDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar Eo = p.Eo(td, sigma_);
    const scalar CdRe = max(this->CdRe(Re), Re*8/3*Eo/(Eo + 4));

    return forceSuSp(Zero, mass*0.75*muc*CdRe/(p.rho()*sqr(p.d())));
}


// ************************************************************************* //
