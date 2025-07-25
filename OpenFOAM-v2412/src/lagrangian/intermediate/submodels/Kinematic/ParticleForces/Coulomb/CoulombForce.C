/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "CoulombForce.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::volVectorField& Foam::CoulombForce<CloudType>::getOrReadField
(
    const word& fieldName
) const
{
    auto* ptr = this->mesh().template getObjectPtr<volVectorField>(fieldName);

    if (!ptr)
    {
        ptr = new volVectorField
        (
            IOobject
            (
                fieldName,
                this->mesh().time().timeName(),
                this->mesh().thisDb(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            this->mesh()
        );
        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CoulombForce<CloudType>::CoulombForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    qPtr_
    (
        Function1<scalar>::New("q", this->coeffs(), &mesh)
    ),
    Ename_(this->coeffs().template getOrDefault<word>("E", "E")),
    EInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::CoulombForce<CloudType>::CoulombForce
(
    const CoulombForce& pf
)
:
    ParticleForce<CloudType>(pf),
    qPtr_(pf.qPtr_.clone()),
    Ename_(pf.Ename_),
    EInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CoulombForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& E = getOrReadField(Ename_);

        EInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                E
            ).ptr()
        );
    }
    else
    {
        EInterpPtr_.reset(nullptr);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::CoulombForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const interpolation<vector>& EInterp = *EInterpPtr_;

    const scalar q = qPtr_->value(p.d());

    // (YSSD:Eq. 6 - left term)
    return forceSuSp
    (
        q*EInterp.interpolate(p.coordinates(), p.currentTetIndices()),
        Zero
    );
}


// ************************************************************************* //
