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

#include "InjectedParticleInjection.H"
#include "mathematicalConstants.H"
#include "bitSet.H"
#include "injectedParticleCloud.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectedParticleInjection<CloudType>::initialise()
{
    const injectedParticleCloud cloud(this->owner().mesh(), cloudName_);

    label nParticles = cloud.size();
    List<scalar> time(nParticles);
    List<vector> position(nParticles);
    List<scalar> diameter(nParticles);
    List<vector> U(nParticles);

    label particlei = 0;

    for (const injectedParticle& p : cloud)
    {
        time[particlei] = p.soi();
        position[particlei] = p.position() + positionOffset_;
        diameter[particlei] = p.d();
        U[particlei] = p.U();

        ++particlei;
    }

    // Combine all proc data
    if (Pstream::parRun())
    {
        List<List<scalar>> procTime(Pstream::nProcs());
        procTime[Pstream::myProcNo()].transfer(time);
        Pstream::allGatherList(procTime);
        time =
            ListListOps::combine<List<scalar>>
            (
                procTime, accessOp<List<scalar>>()
            );

        List<List<point>> procPosition(Pstream::nProcs());
        procPosition[Pstream::myProcNo()].transfer(position);
        Pstream::allGatherList(procPosition);
        position =
            ListListOps::combine<List<point>>
            (
                procPosition, accessOp<List<point>>()
            );

        List<List<scalar>> procD(Pstream::nProcs());
        procD[Pstream::myProcNo()].transfer(diameter);
        Pstream::allGatherList(procD);
        diameter =
            ListListOps::combine<List<scalar>>
            (
                procD, accessOp<List<scalar>>()
            );

        List<List<vector>> procU(Pstream::nProcs());
        procU[Pstream::myProcNo()].transfer(U);
        Pstream::allGatherList(procU);
        U =
            ListListOps::combine<List<vector>>
            (
                procU, accessOp<List<vector>>()
            );
    }

    nParticles = time.size();

    // Reset SOI according to user selection
    scalar minTime = min(time);
    forAll(time, i)
    {
        time[i] -= minTime;
    }

    // Sort and renumber to ensure lists in ascending time
    labelList sortedIndices(Foam::sortedOrder(time));

    time_ = UIndirectList<scalar>(time, sortedIndices);
    position_ = UIndirectList<point>(position, sortedIndices);
    diameter_ = UIndirectList<scalar>(diameter, sortedIndices);
    U_ = UIndirectList<vector>(U, sortedIndices);

    // Pre-calculate injected particle volumes
    List<scalar> volume(nParticles);
    scalar sumVolume = 0;
    forAll(volume, i)
    {
        scalar vol = pow3(diameter_[i])*mathematical::pi/6.0;
        volume[i] = vol;
        sumVolume += vol;
    }
    volume_.transfer(volume);

    // Set the volume of particles to inject
    this->volumeTotal_ = sumVolume;

    // Provide some feedback
    Info<< "    Read " << nParticles << " particles" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectedParticleInjection<CloudType>::InjectedParticleInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    cloudName_(this->coeffDict().lookup("cloud")),
    injectorCells_(),
    injectorTetFaces_(),
    injectorTetPts_(),
    time_(this->template getModelProperty<scalarList>("time")),
    position_(this->template getModelProperty<vectorList>("position")),
    positionOffset_(this->coeffDict().lookup("positionOffset")),
    diameter_(this->template getModelProperty<scalarList>("diameter")),
    U_(this->template getModelProperty<vectorList>("U")),
    volume_(this->template getModelProperty<scalarList>("volume")),
    ignoreOutOfBounds_
    (
        this->coeffDict().getOrDefault("ignoreOutOfBounds", false)
    ),
    currentParticlei_
    (
        this->template getModelProperty<label>
        (
            "currentParticlei",
            -1
        )
    )
{
    if (this->parcelBasis_ != InjectionModel<CloudType>::pbFixed)
    {
        FatalErrorInFunction
            << "Injector model: " << this->modelName()
            << " Parcel basis must be set to fixed"
            << exit(FatalError);
    }

    if (!time_.size())
    {
        // Clean start
        initialise();
    }

    injectorCells_.setSize(position_.size());
    injectorTetFaces_.setSize(position_.size());
    injectorTetPts_.setSize(position_.size());

    updateMesh();

    this->massTotal_ = this->volumeTotal_*this->owner().constProps().rho0();
}


template<class CloudType>
Foam::InjectedParticleInjection<CloudType>::InjectedParticleInjection
(
    const InjectedParticleInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    cloudName_(im.cloudName_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    time_(im.time_),
    position_(im.position_),
    positionOffset_(im.positionOffset_),
    diameter_(im.diameter_),
    U_(im.U_),
    volume_(im.volume_),
    ignoreOutOfBounds_(im.ignoreOutOfBounds_),
    currentParticlei_(im.currentParticlei_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectedParticleInjection<CloudType>::updateMesh()
{
    label nRejected = 0;

    bitSet keep(position_.size(), true);

    forAll(position_, particlei)
    {
        if
        (
            !this->findCellAtPosition
            (
                injectorCells_[particlei],
                injectorTetFaces_[particlei],
                injectorTetPts_[particlei],
                position_[particlei],
                !ignoreOutOfBounds_
            )
        )
        {
            keep.unset(particlei);
            nRejected++;
        }
    }

    if (nRejected > 0)
    {
        inplaceSubset(keep, time_);
        inplaceSubset(keep, position_);
        inplaceSubset(keep, diameter_);
        inplaceSubset(keep, U_);
        inplaceSubset(keep, volume_);
        inplaceSubset(keep, injectorCells_);
        inplaceSubset(keep, injectorTetFaces_);
        inplaceSubset(keep, injectorTetPts_);

        Info<< "    " << nRejected
            << " particles ignored, out of bounds" << endl;
    }
}


template<class CloudType>
Foam::scalar Foam::InjectedParticleInjection<CloudType>::timeEnd() const
{
    return max(time_);
}


template<class CloudType>
Foam::label Foam::InjectedParticleInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    label nParticles = 0;
    forAll(time_, particlei)
    {
        if ((time_[particlei] >= time0) && (time_[particlei] < time1))
        {
            ++nParticles;
        }
    }

    return nParticles;
}


template<class CloudType>
Foam::scalar Foam::InjectedParticleInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar sumVolume = 0;
    forAll(time_, particlei)
    {
        if ((time_[particlei] >= time0) && (time_[particlei] < time1))
        {
            sumVolume += volume_[particlei];
        }
    }

    return sumVolume;
}


template<class CloudType>
void Foam::InjectedParticleInjection<CloudType>::setPositionAndCell
(
    const label parceli,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    // Note: optimisation - consume particle from lists to reduce storage
    // as injection proceeds

    currentParticlei_++;

    position = position_[currentParticlei_];
    cellOwner = injectorCells_[currentParticlei_];
    tetFacei = injectorTetFaces_[currentParticlei_];
    tetPti = injectorTetPts_[currentParticlei_];
}


template<class CloudType>
void Foam::InjectedParticleInjection<CloudType>::setProperties
(
    const label parceli,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity
    parcel.U() = U_[currentParticlei_];

    // Set particle diameter
    parcel.d() = diameter_[currentParticlei_];
}


template<class CloudType>
bool Foam::InjectedParticleInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::InjectedParticleInjection<CloudType>::validInjection
(
    const label
)
{
    return true;
}


template<class CloudType>
void Foam::InjectedParticleInjection<CloudType>::info()
{
    InjectionModel<CloudType>::info();

    if (this->writeTime())
    {
        this->setModelProperty("currentParticlei", currentParticlei_);
        this->setModelProperty("time", time_);
        this->setModelProperty("position", position_);
        this->setModelProperty("diameter", diameter_);
        this->setModelProperty("U", U_);
        this->setModelProperty("volume", volume_);
    }
}


// ************************************************************************* //
