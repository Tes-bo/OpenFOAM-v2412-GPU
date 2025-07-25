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

#include "ConeInjection.H"
#include "Function1.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::ConeInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    positionAxis_(this->coeffDict().lookup("positionAxis")),
    injectorCells_(positionAxis_.size()),
    injectorTetFaces_(positionAxis_.size()),
    injectorTetPts_(positionAxis_.size()),
    duration_(this->coeffDict().getScalar("duration")),
    parcelsPerInjector_
    (
        this->coeffDict().getScalar("parcelsPerInjector")
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
    Umag_
    (
        Function1<scalar>::New
        (
            "Umag",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    thetaInner_
    (
        Function1<scalar>::New
        (
            "thetaInner",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    thetaOuter_
    (
        Function1<scalar>::New
        (
            "thetaOuter",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"), owner.rndGen()
        )
    ),
    nInjected_(Pstream::master() ? this->parcelsAddedTotal() : 0),
    injectorOrder_(identity(positionAxis_.size())),
    tanVec1_(),
    tanVec2_()
{
    updateMesh();

    tanVec1_.setSize(positionAxis_.size());
    tanVec2_.setSize(positionAxis_.size());

    // Convert from user time to reduce the number of time conversion calls
    const Time& time = owner.db().time();
    duration_ = time.userTimeToTime(duration_);
    flowRateProfile_->userTimeToTime(time);
    Umag_->userTimeToTime(time);
    thetaInner_->userTimeToTime(time);
    thetaOuter_->userTimeToTime(time);

    // Normalise direction vector and determine direction vectors
    // tangential to injector axis direction
    forAll(positionAxis_, i)
    {
        vector& axis = positionAxis_[i].second();
        axis.normalise();

        vector tangent = Zero;
        scalar magTangent = 0.0;

        Random& rnd = this->owner().rndGen();
        while (magTangent < SMALL)
        {
            vector v = rnd.sample01<vector>();

            tangent = v - (v & axis)*axis;
            magTangent = mag(tangent);
        }

        tanVec1_[i] = tangent/magTangent;
        tanVec2_[i] = axis^tanVec1_[i];
    }

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_->integrate(0.0, duration_);
}


template<class CloudType>
Foam::ConeInjection<CloudType>::ConeInjection
(
    const ConeInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    positionAxis_(im.positionAxis_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    duration_(im.duration_),
    parcelsPerInjector_(im.parcelsPerInjector_),
    flowRateProfile_(im.flowRateProfile_.clone()),
    Umag_(im.Umag_.clone()),
    thetaInner_(im.thetaInner_.clone()),
    thetaOuter_(im.thetaOuter_.clone()),
    sizeDistribution_(im.sizeDistribution_.clone()),
    nInjected_(im.nInjected_),
    injectorOrder_(im.injectorOrder_),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeInjection<CloudType>::updateMesh()
{
    bitSet reject(positionAxis_.size());

    // Set/cache the injector cells
    forAll(positionAxis_, i)
    {
        if
        (
            !this->findCellAtPosition
            (
                injectorCells_[i],
                injectorTetFaces_[i],
                injectorTetPts_[i],
                positionAxis_[i].first(),
                !this->ignoreOutOfBounds_

            )
        )
        {
            reject.set(i);
        }
    }

    const label nRejected = reject.count();

    if (nRejected)
    {
        reject.flip();
        inplaceSubset(reject, injectorCells_);
        inplaceSubset(reject, injectorTetFaces_);
        inplaceSubset(reject, injectorTetPts_);
        inplaceSubset(reject, positionAxis_);

        Info<< "    " << nRejected
            << " positions rejected, out of bounds" << endl;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::ConeInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        const scalar targetVolume = flowRateProfile_->integrate(0, time1);

        const scalar volumeFraction = targetVolume/this->volumeTotal_;

        const label targetParcels =
            ceil(positionAxis_.size()*parcelsPerInjector_*volumeFraction);

        return targetParcels - returnReduce(nInjected_, sumOp<label>());
    }

    return 0;
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::volumeToInject
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
void Foam::ConeInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    Random& rnd = this->owner().rndGen();
    rnd.shuffle(injectorOrder_);

    const label i = injectorOrder_[parcelI % positionAxis_.size()];
    position = positionAxis_[i].first();
    cellOwner = injectorCells_[i];
    tetFacei = injectorTetFaces_[i];
    tetPti = injectorTetPts_[i];
}


template<class CloudType>
void Foam::ConeInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    Random& rnd = this->owner().rndGen();

    const label i = injectorOrder_[parcelI % positionAxis_.size()];

    // Set direction vectors for position i
    scalar t = time - this->SOI_;
    scalar ti = thetaInner_->value(t);
    scalar to = thetaOuter_->value(t);
    scalar coneAngle = degToRad(rnd.position<scalar>(ti, to));

    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);
    scalar beta = twoPi*rnd.sample01<scalar>();

    vector normal = alpha*(tanVec1_[i]*cos(beta) + tanVec2_[i]*sin(beta));
    vector dirVec = dcorr*positionAxis_[i].second();
    dirVec += normal;
    dirVec.normalise();

    // Set particle velocity
    parcel.U() = Umag_->value(t)*dirVec;

    // Set particle diameter
    parcel.d() = sizeDistribution_().sample();

    // Increment number of particles injected
    nInjected_++;
}


template<class CloudType>
bool Foam::ConeInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ConeInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
