/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "NonInertialFrameForce.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonInertialFrameForce<CloudType>::NonInertialFrameForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    WName_
    (
        this->coeffs().template getOrDefault<word>
        (
            "linearAcceleration",
            "linearAcceleration"
        )
    ),
    W_(Zero),
    omegaName_
    (
        this->coeffs().template getOrDefault<word>
        (
            "angularVelocity",
            "angularVelocity"
        )
    ),
    omega_(Zero),
    omegaDotName_
    (
        this->coeffs().template getOrDefault<word>
        (
            "angularAcceleration",
            "angularAcceleration"
        )
    ),
    omegaDot_(Zero),
    centreOfRotationName_
    (
        this->coeffs().template getOrDefault<word>
        (
            "centreOfRotation",
            "centreOfRotation"
        )
    ),
    centreOfRotation_(Zero)
{}


template<class CloudType>
Foam::NonInertialFrameForce<CloudType>::NonInertialFrameForce
(
    const NonInertialFrameForce& niff
)
:
    ParticleForce<CloudType>(niff),
    WName_(niff.WName_),
    W_(niff.W_),
    omegaName_(niff.omegaName_),
    omega_(niff.omega_),
    omegaDotName_(niff.omegaDotName_),
    omegaDot_(niff.omegaDot_),
    centreOfRotationName_(niff.centreOfRotationName_),
    centreOfRotation_(niff.centreOfRotation_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonInertialFrameForce<CloudType>::~NonInertialFrameForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NonInertialFrameForce<CloudType>::cacheFields(const bool store)
{
    W_ = Zero;
    omega_ = Zero;
    omegaDot_ = Zero;
    centreOfRotation_ = Zero;

    if (store)
    {
        if
        (
            this->mesh().template foundObject<uniformDimensionedVectorField>
            (
                WName_
            )
        )
        {
            const uniformDimensionedVectorField& W = this->mesh().template
                lookupObject<uniformDimensionedVectorField>(WName_);

            W_ = W.value();
        }

        if
        (
            this->mesh().template foundObject<uniformDimensionedVectorField>
            (
                omegaName_
            )
        )
        {
            const uniformDimensionedVectorField& omega = this->mesh().template
                lookupObject<uniformDimensionedVectorField>(omegaName_);

            omega_ = omega.value();
        }

        if
        (
            this->mesh().template foundObject<uniformDimensionedVectorField>
            (
                omegaDotName_
            )
        )
        {
            const uniformDimensionedVectorField& omegaDot =
                this->mesh().template
                lookupObject<uniformDimensionedVectorField>(omegaDotName_);

            omegaDot_ = omegaDot.value();
        }

        if
        (
            this->mesh().template foundObject<uniformDimensionedVectorField>
            (
                centreOfRotationName_
            )
        )
        {
            const uniformDimensionedVectorField& centreOfRotation =
                this->mesh().template
                lookupObject<uniformDimensionedVectorField>
                (
                    centreOfRotationName_
                );

            centreOfRotation_ = centreOfRotation.value();
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::NonInertialFrameForce<CloudType>::calcNonCoupled
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

    const vector r = p.position() - centreOfRotation_;

    value.Su() =
        mass
       *(
           -W_
          + (r ^ omegaDot_)
          + 2.0*(p.U() ^ omega_)
          + (omega_ ^ (r ^ omega_))
        );

    return value;
}


// ************************************************************************* //
