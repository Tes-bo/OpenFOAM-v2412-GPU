/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "PatchFlowRateInjection.H"
#include "distributionModel.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchFlowRateInjection<CloudType>::PatchFlowRateInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName,typeName),
    patchInjectionBase(owner.mesh(), this->coeffDict().getWord("patch")),
    phiName_(this->coeffDict().template getOrDefault<word>("phi", "phi")),
    rhoName_(this->coeffDict().template getOrDefault<word>("rho", "rho")),
    duration_(this->coeffDict().getScalar("duration")),
    concentration_
    (
        Function1<scalar>::New
        (
            "concentration",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    parcelConcentration_
    (
        this->coeffDict().getScalar("parcelConcentration")
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    )
{
    // Convert from user time to reduce the number of time conversion calls
    const Time& time = owner.db().time();
    duration_ = time.userTimeToTime(duration_);
    concentration_->userTimeToTime(time);

    patchInjectionBase::updateMesh(owner.mesh());

    // Re-initialise total mass/volume to inject to zero
    // - will be reset during each injection
    this->volumeTotal_ = 0.0;
    this->massTotal_ = 0.0;
}


template<class CloudType>
Foam::PatchFlowRateInjection<CloudType>::PatchFlowRateInjection
(
    const PatchFlowRateInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    patchInjectionBase(im),
    phiName_(im.phiName_),
    rhoName_(im.rhoName_),
    duration_(im.duration_),
    concentration_(im.concentration_.clone()),
    parcelConcentration_(im.parcelConcentration_),
    sizeDistribution_(im.sizeDistribution_.clone())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchFlowRateInjection<CloudType>::~PatchFlowRateInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchFlowRateInjection<CloudType>::updateMesh()
{
    patchInjectionBase::updateMesh(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::flowRate() const
{
   const polyMesh& mesh = this->owner().mesh();

    const auto& phi = mesh.lookupObject<surfaceScalarField>(phiName_);

    const scalarField& phip = phi.boundaryField()[patchId_];

    scalar flowRateIn = 0.0;
    if (phi.dimensions() == dimVolume/dimTime)
    {
        flowRateIn = max(0.0, -sum(phip));
    }
    else
    {
        const auto& rho = mesh.lookupObject<volScalarField>(rhoName_);
        const scalarField& rhop = rho.boundaryField()[patchId_];

        flowRateIn = max(0.0, -sum(phip/rhop));
    }

    reduce(flowRateIn, sumOp<scalar>());

    return flowRateIn;
}


template<class CloudType>
Foam::label Foam::PatchFlowRateInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar dt = time1 - time0;

        scalar c = concentration_->value(0.5*(time0 + time1));

        scalar nParcels = parcelConcentration_*c*flowRate()*dt;

        Random& rnd = this->owner().rndGen();

        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        if
        (
            nParcelsToInject > 0
         && (
               nParcels - scalar(nParcelsToInject)
             > rnd.globalPosition(scalar(0), scalar(1))
            )
        )
        {
            ++nParcelsToInject;
        }

        return nParcelsToInject;
    }

    return 0;
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar volume = 0.0;

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar c = concentration_->value(0.5*(time0 + time1));

        volume = c*(time1 - time0)*flowRate();
    }

    this->volumeTotal_ = volume;
    this->massTotal_ = volume*this->owner().constProps().rho0();

    return volume;
}


template<class CloudType>
void Foam::PatchFlowRateInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    patchInjectionBase::setPositionAndCell
    (
        this->owner().mesh(),
        this->owner().rndGen(),
        position,
        cellOwner,
        tetFacei,
        tetPti
    );
}


template<class CloudType>
void Foam::PatchFlowRateInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity to carrier velocity
    parcel.U() = this->owner().U()[parcel.cell()];

    // Set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::PatchFlowRateInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::PatchFlowRateInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
