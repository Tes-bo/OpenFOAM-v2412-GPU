/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "InjectionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModelList<CloudType>::InjectionModelList(CloudType& owner)
:
    PtrList<InjectionModel<CloudType>>()
{}


template<class CloudType>
Foam::InjectionModelList<CloudType>::InjectionModelList
(
    const dictionary& dict,
    CloudType& owner
)
:
    PtrList<InjectionModel<CloudType>>()
{
    Info<< "Constructing particle injection models" << endl;

    label count = dict.size();
    if (count)
    {
        this->resize(count);
    }

    count = 0;
    for (const entry& dEntry : dict)
    {
        const word& model = dEntry.keyword();
        const dictionary& props = dEntry.dict();

        Info<< "Creating injector: " << model << endl;

        this->set
        (
            count,
            InjectionModel<CloudType>::New
            (
                props,
                model,
                props.get<word>("type"),
                owner
            )
        );

        ++count;
    }

    if (!count)
    {
        this->resize(1);

        this->set
        (
            0,
            InjectionModel<CloudType>::New
            (
                dict,
                "none",
                "none",
                owner
            )
        );
    }
}


template<class CloudType>
Foam::InjectionModelList<CloudType>::InjectionModelList
(
    const InjectionModelList<CloudType>& iml
)
:
    PtrList<InjectionModel<CloudType>>(iml)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::InjectionModelList<CloudType>::timeStart() const
{
    scalar minTime = GREAT;
    for (const auto& model : *this)
    {
        minTime = min(minTime, model.timeStart());
    }

    return minTime;
}


template<class CloudType>
Foam::scalar Foam::InjectionModelList<CloudType>::timeEnd() const
{
    scalar maxTime = -GREAT;
    for (const auto& model : *this)
    {
        maxTime = max(maxTime, model.timeEnd());
    }

    return maxTime;
}


template<class CloudType>
Foam::scalar Foam::InjectionModelList<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar vol = 0.0;
    for (auto& model : *this)
    {
        vol += model.volumeToInject(time0, time1);
    }

    return vol;
}


template<class CloudType>
Foam::scalar Foam::InjectionModelList<CloudType>::averageParcelMass()
{
    scalar mass = 0.0;
    scalar massTotal = 0.0;
    for (auto& model : *this)
    {
        scalar mt = model.massTotal();
        mass += mt*model.averageParcelMass();
        massTotal += mt;
    }

    return mass/massTotal;
}


template<class CloudType>
void Foam::InjectionModelList<CloudType>::updateMesh()
{
    for (auto& model : *this)
    {
        model.updateMesh();
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::InjectionModelList<CloudType>::inject
(
    TrackCloudType& cloud,
    typename CloudType::parcelType::trackingData& td
)
{
    for (auto& model : *this)
    {
        model.inject(cloud, td);
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::InjectionModelList<CloudType>::injectSteadyState
(
    TrackCloudType& cloud,
    typename CloudType::parcelType::trackingData& td,
    const scalar trackTime
)
{
    for (auto& model : *this)
    {
        model.injectSteadyState(cloud, td, trackTime);
    }
}


template<class CloudType>
void Foam::InjectionModelList<CloudType>::info()
{
    for (auto& model : *this)
    {
        model.info();
    }
}


// ************************************************************************* //
