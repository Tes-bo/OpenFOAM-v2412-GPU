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

#include "VoidFraction.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFraction<CloudType>::write()
{
    if (thetaPtr_)
    {
        thetaPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "thetaPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFraction<CloudType>::VoidFraction
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    thetaPtr_(nullptr)
{}


template<class CloudType>
Foam::VoidFraction<CloudType>::VoidFraction
(
    const VoidFraction<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    thetaPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFraction<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    if (thetaPtr_)
    {
        thetaPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        thetaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Theta",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless, Zero)
            )
        );
    }
}


template<class CloudType>
void Foam::VoidFraction<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    volScalarField& theta = thetaPtr_();

    const fvMesh& mesh = this->owner().mesh();

    theta.primitiveFieldRef() /= mesh.time().deltaTValue()*mesh.V();

    CloudFunctionObject<CloudType>::postEvolve(td);
}


template<class CloudType>
bool Foam::VoidFraction<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point&,
    const typename parcelType::trackingData& td
)
{
    volScalarField& theta = thetaPtr_();

    theta[p.cell()] += dt*p.nParticle()*p.volume();

    return true;
}


// ************************************************************************* //
