/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "PatchInjection.H"
#include "distributionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::Enum<typename Foam::PatchInjection<CloudType>::velocityType>
Foam::PatchInjection<CloudType>::velocityTypeNames_
({
    { vtFixedValue, "fixedValue" },
    { vtPatchValue, "patchValue" },
    { vtZeroGradient, "zeroGradient" },
});

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInjection<CloudType>::PatchInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    patchInjectionBase(owner.mesh(), this->coeffDict().getWord("patch")),
    duration_(this->coeffDict().getScalar("duration")),
    parcelsPerSecond_
    (
        this->coeffDict().getScalar("parcelsPerSecond")
    ),
    velocityType_
    (
        velocityTypeNames_.getOrDefault
        (
            "velocityType",
            this->coeffDict(),
            vtFixedValue
        )
    ),
    U0_
    (
        velocityType_ == vtFixedValue
      ? this->coeffDict().template get<vector>("U0")
      : Zero
    ),
    flowRateProfile_
    (
        Function1<scalar>::New
        (
            "flowRateProfile",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    currentParceli_(-1),
    currentFacei_(-1)
{
    // Convert from user time to reduce the number of time conversion calls
    const Time& time = owner.db().time();
    duration_ = time.userTimeToTime(duration_);
    flowRateProfile_->userTimeToTime(time);

    patchInjectionBase::updateMesh(owner.mesh());

    // Set total volume/mass to inject
    this->volumeTotal_ = flowRateProfile_->integrate(0.0, duration_);
}


template<class CloudType>
Foam::PatchInjection<CloudType>::PatchInjection
(
    const PatchInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    patchInjectionBase(im),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    velocityType_(im.velocityType_),
    U0_(im.U0_),
    flowRateProfile_(im.flowRateProfile_.clone()),
    sizeDistribution_(im.sizeDistribution_.clone()),
    currentParceli_(im.currentParceli_),
    currentFacei_(im.currentFacei_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchInjection<CloudType>::updateMesh()
{
    patchInjectionBase::updateMesh(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::PatchInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::PatchInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar nParcels = (time1 - time0)*parcelsPerSecond_;
        Random& rnd = this->owner().rndGen();
        scalar rndPos = rnd.globalPosition(scalar(0), scalar(1));
        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        if
        (
            nParcelsToInject > 0
         && (nParcels - scalar(nParcelsToInject) > rndPos)
        )
        {
            ++nParcelsToInject;
        }

        return nParcelsToInject;
    }

    return 0;
}


template<class CloudType>
Foam::scalar Foam::PatchInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return flowRateProfile_->integrate(time0, time1);
    }

    return 0.0;
}


template<class CloudType>
void Foam::PatchInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    currentParceli_ = parcelI;

    currentFacei_ = patchInjectionBase::setPositionAndCell
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
void Foam::PatchInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity
    switch (velocityType_)
    {
        case vtFixedValue:
        {
            parcel.U() = U0_;
            break;
        }
        case vtPatchValue:
        {
            if (parcelI != currentParceli_)
            {
                WarningInFunction
                    << "Synchronisation problem: "
                    << "attempting to set injected parcel " << parcelI
                    << " properties using cached parcel " << currentParceli_
                    << " properties" << endl;
            }

            const label patchFacei = currentFacei_;

            if (patchFacei < 0)
            {
                FatalErrorInFunction
                    << "Unable to set parcel velocity using patch value "
                    << "due to missing face index: patchFacei=" << patchFacei
                    << abort(FatalError);
            }

            const volVectorField& U = this->owner().U();
            const label patchi = this->patchId_;
            parcel.U() = U.boundaryField()[patchi][patchFacei];
            break;
        }
        case vtZeroGradient:
        {
            const label celli = parcel.cell();

            if (celli < 0)
            {
                FatalErrorInFunction
                    << "Unable to set parcel velocity using zeroGradient "
                    << "due to missing cell index"
                    << abort(FatalError);
            }

            const volVectorField& U = this->owner().U();
            parcel.U() = U[celli];
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled velocityType "
                << velocityTypeNames_[velocityType_]
                << ". Available options are:"
                << velocityTypeNames_.sortedToc()
                << abort(FatalError);
        }
    }

    // Set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::PatchInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::PatchInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
