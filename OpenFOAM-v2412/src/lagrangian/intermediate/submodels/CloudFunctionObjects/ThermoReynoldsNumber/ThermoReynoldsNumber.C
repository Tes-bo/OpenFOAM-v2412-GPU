/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "ThermoReynoldsNumber.H"
#include "ThermoCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoReynoldsNumber<CloudType>::ThermoReynoldsNumber
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName)
{}


template<class CloudType>
Foam::ThermoReynoldsNumber<CloudType>::ThermoReynoldsNumber
(
    const ThermoReynoldsNumber<CloudType>& re
)
:
    CloudFunctionObject<CloudType>(re)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoReynoldsNumber<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    auto& c = this->owner();

    auto* resultPtr = c.template getObjectPtr<IOField<scalar>>("Re");

    if (!resultPtr)
    {
        resultPtr = new IOField<scalar>
        (
            IOobject
            (
                "Re",
                c.time().timeName(),
                c,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        );

        resultPtr->store();
    }
    auto& Re = *resultPtr;

    Re.resize(c.size());

    typename parcelType::trackingData& nctd =
        const_cast<typename parcelType::trackingData&>(td);

    label parceli = 0;
    for (const parcelType& p : c)
    {
        scalar Ts, rhos, mus, Pr, kappas;
        p.template calcSurfaceValues<CloudType>
        (
            c, nctd, p.T(), Ts, rhos, mus, Pr, kappas
        );

        Re[parceli++] = p.Re(rhos, p.U(), td.Uc(), p.d(), mus);
    }

    const bool haveParticles = c.size();
    if (c.time().writeTime() && returnReduceOr(haveParticles))
    {
        Re.write(haveParticles);
    }
}


// ************************************************************************* //
