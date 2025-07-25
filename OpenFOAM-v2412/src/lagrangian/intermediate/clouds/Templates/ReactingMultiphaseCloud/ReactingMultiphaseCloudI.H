/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ReactingMultiphaseCloud<CloudType>&
Foam::ReactingMultiphaseCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ReactingMultiphaseCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ReactingMultiphaseCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::DevolatilisationModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::devolatilisation() const
{
    return devolatilisationModel_;
}


template<class CloudType>
inline Foam::DevolatilisationModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::devolatilisation()
{
    return *devolatilisationModel_;
}


template<class CloudType>
inline const Foam::SurfaceReactionModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::surfaceReaction() const
{
    return surfaceReactionModel_;
}


template<class CloudType>
inline Foam::SurfaceReactionModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::surfaceReaction()
{
    return *surfaceReactionModel_;
}


template<class CloudType>
inline void Foam::ReactingMultiphaseCloud<CloudType>::transferToCarrier
(
    const parcelType& p,
    const typename parcelType::trackingData& td
)
{
    const label celli = p.cell();

    const scalar m = p.nParticle()*p.mass();
    const scalar pc = td.pc();
    const scalar T = p.T();

    this->rhokTrans()[celli] += m;

    this->UTrans()[celli] += m*p.U();

    const auto& comp = this->composition();
    const label idG = comp.idGas();
    const label idL = comp.idLiquid();
    // const label idS = composition.idSolid();

    // Absorb parcel into carrier phase
    auto transfer = [&]
    (
        const label phaseType,
        const label phasei,
        const scalarField& Y
    )
    {
        const scalar YMix = p.Y()[phaseType];

        forAll(Y, i)
        {
            const label gid = comp.localToCarrierId(phaseType, i);
            this->rhoTrans(gid)[celli] += m*YMix*Y[i];
            this->hsTrans()[celli] += m*YMix*comp.Hs(phasei, Y, pc, T);
        }
    };

    transfer(parcelType::GAS, idG, p.YGas());
    transfer(parcelType::LIQ, idL, p.YLiquid());
    // No mapping between solid components and carrier phase
    //transfer(parcelType::SLD, idS, p.YSolid());
}


// ************************************************************************* //
