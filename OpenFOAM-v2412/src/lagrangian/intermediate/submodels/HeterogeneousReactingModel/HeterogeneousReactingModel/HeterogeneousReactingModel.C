/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "HeterogeneousReactingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HeterogeneousReactingModel<CloudType>::HeterogeneousReactingModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner),
    dMass_(0.0),
    nF_(0)
{}


template<class CloudType>
Foam::HeterogeneousReactingModel<CloudType>::HeterogeneousReactingModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    dMass_(0.0),
    nF_(this->coeffDict().template getOrDefault<label>("nF", 1))
{}


template<class CloudType>
Foam::HeterogeneousReactingModel<CloudType>::HeterogeneousReactingModel
(
    const HeterogeneousReactingModel<CloudType>& srm
)
:
    CloudSubModelBase<CloudType>(srm),
    dMass_(srm.dMass_),
    nF_(srm.nF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::HeterogeneousReactingModel<CloudType>::addToSurfaceReactionMass
(
    const scalar dMass
)
{
    dMass_ += dMass;
}


template<class CloudType>
Foam::label Foam::HeterogeneousReactingModel<CloudType>::nF() const
{
    return nF_;
}


template<class CloudType>
void Foam::HeterogeneousReactingModel<CloudType>::info()
{
    CloudSubModelBase<CloudType>::info();

    const scalar mass0 = this->template getBaseProperty<scalar>("mass");
    const scalar massTotal = mass0 + returnReduce(dMass_, sumOp<scalar>());

    Log_<< "    Mass transfer surface reaction  = " << massTotal << nl;

    if (this->writeTime())
    {
        this->setBaseProperty("mass", massTotal);
        dMass_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HeterogeneousReactingModelNew.C"

// ************************************************************************* //
