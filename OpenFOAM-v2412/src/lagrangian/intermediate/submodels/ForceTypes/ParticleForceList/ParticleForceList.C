/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "ParticleForceList.H"
#include "entry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleForceList<CloudType>::ParticleForceList
(
    CloudType& owner,
    const fvMesh& mesh
)
:
    PtrList<ParticleForce<CloudType>>(),
    owner_(owner),
    mesh_(mesh),
    dict_(),
    calcCoupled_(true),
    calcNonCoupled_(true)
{}


template<class CloudType>
Foam::ParticleForceList<CloudType>::ParticleForceList
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool readFields
)
:
    PtrList<ParticleForce<CloudType>>(),
    owner_(owner),
    mesh_(mesh),
    dict_(dict),
    calcCoupled_(true),
    calcNonCoupled_(true)
{
    if (readFields)
    {
        Info<< "Constructing particle forces" << endl;

        this->resize(dict.size());

        label count = 0;
        for (const entry& dEntry : dict)
        {
            const word& modelName = dEntry.keyword();
            if (dEntry.isDict())
            {
                this->set
                (
                    count,
                    ParticleForce<CloudType>::New
                    (
                        owner,
                        mesh,
                        dEntry.dict(),
                        modelName
                    )
                );
            }
            else
            {
                this->set
                (
                    count,
                    ParticleForce<CloudType>::New
                    (
                        owner,
                        mesh,
                        dict,
                        modelName
                    )
                );
            }
            ++count;
        }

        if (!count)
        {
            Info<< "    none" << endl;
        }
    }
}


template<class CloudType>
Foam::ParticleForceList<CloudType>::ParticleForceList
(
    const ParticleForceList& pf
)
:
    PtrList<ParticleForce<CloudType>>(pf),
    owner_(pf.owner_),
    mesh_(pf.mesh_),
    dict_(pf.dict_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleForceList<CloudType>::~ParticleForceList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleForceList<CloudType>::cacheFields(const bool store)
{
    forAll(*this, i)
    {
        this->operator[](i).cacheFields(store);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::ParticleForceList<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero);

    if (calcCoupled_)
    {
        forAll(*this, i)
        {
            value += this->operator[](i).calcCoupled(p, td, dt, mass, Re, muc);
        }
    }

    return value;
}


template<class CloudType>
Foam::forceSuSp Foam::ParticleForceList<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero);

    if (calcNonCoupled_)
    {
        forAll(*this, i)
        {
            value +=
                this->operator[](i).calcNonCoupled(p, td, dt, mass, Re, muc);
        }
    }

    return value;
}


template<class CloudType>
Foam::scalar Foam::ParticleForceList<CloudType>::massEff
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar mass
) const
{
    scalar massEff = mass;
    forAll(*this, i)
    {
        massEff += this->operator[](i).massAdd(p, td, mass);
    }

    return massEff;
}


// ************************************************************************* //
