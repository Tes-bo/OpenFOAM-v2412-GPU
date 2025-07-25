/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "Rebound.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Rebound<CloudType>::Rebound
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    UFactor_(this->coeffDict().getScalar("UFactor"))
{}


template<class CloudType>
Foam::Rebound<CloudType>::Rebound(const Rebound<CloudType>& pim)
:
    PatchInteractionModel<CloudType>(pim),
    UFactor_(pim.UFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::Rebound<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    vector& U = p.U();

    keepParticle = true;
    p.active(true);

    vector nw;
    vector Up;

    this->owner().patchData(p, pp, nw, Up);

    // Calculate motion relative to patch velocity
    U -= Up;

    scalar Un = U & nw;

    if (Un > 0.0)
    {
        U -= UFactor_*2.0*Un*nw;
    }

    // Return velocity to global space
    U += Up;

    return true;
}


// ************************************************************************* //
