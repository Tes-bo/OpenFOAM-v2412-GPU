/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "KinematicWeberNumber.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KinematicWeberNumber<CloudType>::KinematicWeberNumber
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    sigma_(dict.getScalar("sigma"))
{}


template<class CloudType>
Foam::KinematicWeberNumber<CloudType>::KinematicWeberNumber
(
    const KinematicWeberNumber<CloudType>& we
)
:
    CloudFunctionObject<CloudType>(we),
    sigma_(we.sigma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::KinematicWeberNumber<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    auto& c = this->owner();

    auto* resultPtr = c.template getObjectPtr<IOField<scalar>>("We");

    if (!resultPtr)
    {
        resultPtr = new IOField<scalar>
        (
            IOobject
            (
                "We",
                c.time().timeName(),
                c,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        );

        resultPtr->store();
    }
    auto& We = *resultPtr;

    We.resize(c.size());

    label parceli = 0;
    for (const parcelType& p : c)
    {
        We[parceli++] = p.We(td, sigma_);
    }

    const bool haveParticles = c.size();
    if (c.time().writeTime() && returnReduceOr(haveParticles))
    {
        We.write(haveParticles);
    }
}


// ************************************************************************* //
