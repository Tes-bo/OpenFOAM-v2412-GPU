/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "InjectedParticleDistributionInjection.H"
#include "mathematicalConstants.H"
#include "bitSet.H"
#include "injectedParticleCloud.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectedParticleDistributionInjection<CloudType>::initialise()
{
    injectedParticleCloud ipCloud(this->owner().mesh(), cloudName_);

    List<label> tag(ipCloud.size());
    List<point> position(ipCloud.size());
    List<vector> U(ipCloud.size());
    List<scalar> soi(ipCloud.size());
    List<scalar> d(ipCloud.size());

    // Flatten all data
    label particlei = 0;
    for (const injectedParticle& p : ipCloud)
    {
        tag[particlei] = p.tag();
        position[particlei] = p.position();
        U[particlei] = p.U();
        soi[particlei] = p.soi();
        d[particlei] = p.d();
        particlei++;
    }

    // Combine all proc data
    if (Pstream::parRun())
    {
        List<List<label>> procTag(Pstream::nProcs());
        procTag[Pstream::myProcNo()].transfer(tag);
        Pstream::allGatherList(procTag);
        tag =
            ListListOps::combine<List<label>>
            (
                procTag, accessOp<List<label>>()
            );

        List<List<point>> procPosition(Pstream::nProcs());
        procPosition[Pstream::myProcNo()].transfer(position);
        Pstream::allGatherList(procPosition);
        position =
            ListListOps::combine<List<point>>
            (
                procPosition, accessOp<List<point>>()
            );

        List<List<vector>> procU(Pstream::nProcs());
        procU[Pstream::myProcNo()].transfer(U);
        Pstream::allGatherList(procU);
        U =
            ListListOps::combine<List<vector>>
            (
                procU, accessOp<List<vector>>()
            );

        List<List<scalar>> procSOI(Pstream::nProcs());
        procSOI[Pstream::myProcNo()].transfer(soi);
        Pstream::allGatherList(procSOI);
        soi =
            ListListOps::combine<List<scalar>>
            (
                procSOI, accessOp<List<scalar>>()
            );

        List<List<scalar>> procD(Pstream::nProcs());
        procD[Pstream::myProcNo()].transfer(d);
        Pstream::allGatherList(procD);
        d =
            ListListOps::combine<List<scalar>>
            (
                procD, accessOp<List<scalar>>()
            );
    }

    label maxTag = -1;
    forAll(tag, particlei)
    {
        maxTag = max(maxTag, tag[particlei]);
    }

    label nInjectors = maxTag + 1;
    List<scalar> injStartTime(nInjectors, GREAT);
    List<scalar> injEndTime(nInjectors, -GREAT);
    List<DynamicList<point>> injPosition(nInjectors);
    List<DynamicList<vector>> injU(nInjectors);
    List<DynamicList<scalar>> injDiameter(nInjectors);

    // Cache the particle information per tag
    forAll(tag, i)
    {
        const label tagi = tag[i];
        const scalar t = soi[i];
        injStartTime[tagi] = min(t, injStartTime[tagi]);
        injEndTime[tagi] = max(t, injEndTime[tagi]);
        injPosition[tagi].append(position[i]);
        injU[tagi].append(U[i]);
        injDiameter[tagi].append(d[i]);
    }

    // Remove single particles and injectors where injection interval is 0
    // - cannot generate a volume flow rate
    scalar sumVolume = 0;
    startTime_.setSize(nInjectors, 0);
    endTime_.setSize(nInjectors, 0);
    sizeDistribution_.setSize(nInjectors);
    position_.setSize(nInjectors);
    U_.setSize(nInjectors);
    volumeFlowRate_.setSize(nInjectors, 0);

    scalar minTime = GREAT;

    // Populate injector properties, filtering out invalid entries
    Random& rnd = this->owner().rndGen();
    label injectori = 0;
    forAll(injDiameter, i)
    {
        const DynamicList<scalar>& diameters = injDiameter[i];
        const label nParticle = diameters.size();
        const scalar dTime = injEndTime[i] - injStartTime[i];

        if ((nParticle > 1) && (dTime > ROOTVSMALL))
        {
            minTime = min(minTime, injStartTime[i]);

            startTime_[injectori] = injStartTime[i];
            endTime_[injectori] = injEndTime[i];

            // Re-sample the cloud data
            position_[injectori].setSize(resampleSize_);
            U_[injectori].setSize(resampleSize_);
            List<point>& positioni = position_[injectori];
            List<vector>& Ui = U_[injectori];

            for (label samplei = 0; samplei < resampleSize_; ++samplei)
            {
                label posi = rnd.globalPosition<label>(0, nParticle - 1);
                positioni[samplei] = injPosition[i][posi] + positionOffset_;
                Ui[samplei] = injU[i][posi];
            }

            // Calculate the volume flow rate
            scalar sumPow3 = 0;
            forAll(diameters, particlei)
            {
                sumPow3 += pow3(diameters[particlei]);
            }

            const scalar volume = sumPow3*mathematical::pi/6.0;
            sumVolume += volume;
            volumeFlowRate_[injectori] = volume/dTime;

            // Create the size distribution using the user-specified bin width
            sizeDistribution_.set
            (
                injectori,
                new distributionModels::general
                (
                    diameters,
                    binWidth_,
                    this->owner().rndGen()
                )
            );

            injectori++;
        }
    }

    // Resize
    startTime_.setSize(injectori);
    endTime_.setSize(injectori);
    position_.setSize(injectori);
    U_.setSize(injectori);
    volumeFlowRate_.setSize(injectori);
    sizeDistribution_.setSize(injectori);

    // Reset start time to zero
    forAll(startTime_, injectori)
    {
        startTime_[injectori] -= minTime;
        endTime_[injectori] -= minTime;
    }


    // Set the volume of parcels to inject
    this->volumeTotal_ = sumVolume;

    // Provide some feedback
    Info<< "    Read " << position_.size() << " injectors with "
        << tag.size() << " total particles" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectedParticleDistributionInjection<CloudType>::
InjectedParticleDistributionInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    cloudName_(this->coeffDict().lookup("cloud")),
    startTime_(this->template getModelProperty<scalarList>("startTime")),
    endTime_(this->template getModelProperty<scalarList>("endTime")),
    position_(this->template getModelProperty<List<vectorList>>("position")),
    positionOffset_(this->coeffDict().lookup("positionOffset")),
    volumeFlowRate_
    (
        this->template getModelProperty<scalarList>("volumeFlowRate")
    ),
    U_(this->template getModelProperty<List<vectorList>>("U")),
    binWidth_(this->coeffDict().getScalar("binWidth")),
    sizeDistribution_(),
    parcelsPerInjector_
    (
        ceil(this->coeffDict().getScalar("parcelsPerInjector"))
    ),
    resampleSize_
    (
        this->coeffDict().getOrDefault("resampleSize", label(100))
    ),
    applyDistributionMassTotal_
    (
        this->coeffDict().getBool("applyDistributionMassTotal")
    ),
    ignoreOutOfBounds_
    (
        this->coeffDict().getOrDefault("ignoreOutOfBounds", false)
    ),
    nParcelsInjected_(this->parcelsAddedTotal()),
    nParcelsInjected0_(0),
    currentInjectori_(0),
    currentSamplei_(0)
{
    if (startTime_.size())
    {
        // Restart
        sizeDistribution_.setSize(startTime_.size());
        forAll(sizeDistribution_, i)
        {
            const word dictName("distribution" + Foam::name(i));
            dictionary dict;
            this->getModelDict(dictName, dict);

            sizeDistribution_.set
            (
                i,
                new distributionModels::general(dict, this->owner().rndGen())
            );
        }
    }
    else
    {
        // Clean start
        initialise();
    }

    // Set the mass of parcels to inject from distribution if requested
    if (applyDistributionMassTotal_)
    {
        this->massTotal_ = this->volumeTotal_*this->owner().constProps().rho0();
        Info<< "    Set mass to inject from distribution: "
            << this->massTotal_ << endl;
    }
}


template<class CloudType>
Foam::InjectedParticleDistributionInjection<CloudType>::
InjectedParticleDistributionInjection
(
    const InjectedParticleDistributionInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    cloudName_(im.cloudName_),
    startTime_(im.startTime_),
    endTime_(im.endTime_),
    position_(im.position_),
    positionOffset_(im.positionOffset_),
    volumeFlowRate_(im.volumeFlowRate_),
    U_(im.U_),
    binWidth_(im.binWidth_),
    sizeDistribution_(im.sizeDistribution_.size()),
    parcelsPerInjector_(im.parcelsPerInjector_),
    resampleSize_(im.resampleSize_),
    applyDistributionMassTotal_(im.applyDistributionMassTotal_),
    ignoreOutOfBounds_(im.ignoreOutOfBounds_),
    nParcelsInjected_(im.nParcelsInjected_),
    nParcelsInjected0_(im.nParcelsInjected0_),
    currentInjectori_(0),
    currentSamplei_(0)
{
    forAll(sizeDistribution_, injectori)
    {
        if (sizeDistribution_.set(injectori))
        {
            sizeDistribution_.set
            (
                injectori,
                new distributionModels::general
                (
                    im.sizeDistribution_[injectori]
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectedParticleDistributionInjection<CloudType>::
~InjectedParticleDistributionInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectedParticleDistributionInjection<CloudType>::updateMesh()
{}


template<class CloudType>
Foam::scalar
Foam::InjectedParticleDistributionInjection<CloudType>::timeEnd() const
{
    return max(endTime_);
}


template<class CloudType>
Foam::label
Foam::InjectedParticleDistributionInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    // Ensure all procs know the latest parcel count
    nParcelsInjected_ += returnReduce(nParcelsInjected0_, sumOp<label>());
    nParcelsInjected0_ = 0;

    if (startTime_.empty() || this->volumeTotal_ < ROOTVSMALL)
    {
        return 0;
    }

    scalar targetVolume = 0;
    forAll(startTime_, injectori)
    {
        if (time1 > startTime_[injectori])
        {
            scalar totalDuration =
                min(time1, endTime_[injectori]) - startTime_[injectori];

            targetVolume += volumeFlowRate_[injectori]*totalDuration;
        }
    }

    const label targetParcels =
        round
        (
            scalar(startTime_.size()*parcelsPerInjector_)
           *targetVolume/this->volumeTotal_
        );

    const label nParcels = targetParcels - nParcelsInjected_;

    return nParcels;
}


template<class CloudType>
Foam::scalar
Foam::InjectedParticleDistributionInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar volume = 0;
    forAll(startTime_, injectori)
    {
        if ((time1 > startTime_[injectori]) && (time1 <= endTime_[injectori]))
        {
            scalar duration = min(time1, endTime_[injectori]) - time0;
            volume += volumeFlowRate_[injectori]*duration;
        }
    }

    return volume;
}


template<class CloudType>
void Foam::InjectedParticleDistributionInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    Random& rnd = this->owner().rndGen();
    currentInjectori_ = rnd.globalPosition<label>(0, position_.size() - 1);
    currentSamplei_ = rnd.globalPosition<label>(0, resampleSize_ - 1);

    position = position_[currentInjectori_][currentSamplei_];

    // Cache all mesh props for each position?
    this->findCellAtPosition
    (
        cellOwner,
        tetFaceI,
        tetPtI,
        position
    );
}


template<class CloudType>
void Foam::InjectedParticleDistributionInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity
    parcel.U() = U_[currentInjectori_][currentSamplei_];

    // Set particle diameter
    parcel.d() = sizeDistribution_[currentInjectori_].sample();

    // Increment number of particles injected
    // Note: local processor only!
    nParcelsInjected0_++;
}


template<class CloudType>
bool
Foam::InjectedParticleDistributionInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::InjectedParticleDistributionInjection<CloudType>::validInjection
(
    const label
)
{
    return true;
}


template<class CloudType>
void Foam::InjectedParticleDistributionInjection<CloudType>::info()
{
    InjectionModel<CloudType>::info();

    if (this->writeTime())
    {
        this->setModelProperty("startTime", startTime_);
        this->setModelProperty("endTime", endTime_);
        this->setModelProperty("position", position_);
        this->setModelProperty("volumeFlowRate", volumeFlowRate_);
        this->setModelProperty("U", U_);
        forAll(sizeDistribution_, i)
        {
            const distributionModels::general& dist = sizeDistribution_[i];
            const word dictName("distribution" + Foam::name(i));
            dictionary dict(dist.writeDict(dictName));
            this->setModelProperty(dictName, dict);
        }
    }
}


// ************************************************************************* //
