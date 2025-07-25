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

#include "NoheterogeneousReacting.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoheterogeneousReacting<CloudType>::NoheterogeneousReacting
(
    const dictionary&,
    CloudType& owner
)
:
    HeterogeneousReactingModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoheterogeneousReacting<CloudType>::NoheterogeneousReacting
(
    const NoheterogeneousReacting<CloudType>& srm
)
:
    HeterogeneousReactingModel<CloudType>(srm.owner_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoheterogeneousReacting<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::NoheterogeneousReacting<CloudType>::calculate
(
    const scalar,
    const scalar,
    const scalar,
    const label,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalarField&,
    scalarField&,
    const scalar,
    scalar&,
    scalarField&,
    scalarField&
) const
{
    return 0;
}


template<class CloudType>
Foam::label Foam::NoheterogeneousReacting<CloudType>::nReactions() const
{
    return 0;
}


// ************************************************************************* //
