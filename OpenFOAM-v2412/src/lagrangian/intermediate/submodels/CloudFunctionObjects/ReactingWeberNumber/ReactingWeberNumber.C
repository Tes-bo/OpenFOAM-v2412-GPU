/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "ReactingWeberNumber.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingWeberNumber<CloudType>::ReactingWeberNumber
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName)
{}


template<class CloudType>
Foam::ReactingWeberNumber<CloudType>::ReactingWeberNumber
(
    const ReactingWeberNumber<CloudType>& we
)
:
    CloudFunctionObject<CloudType>(we)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingWeberNumber<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    const auto& c = this->owner();

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

    const auto& thermo = c.db().template lookupObject<SLGThermo>("SLGThermo");
    const auto& liquids = thermo.liquids();

    const auto& UInterp = td.UInterp();
    const auto& pInterp = td.pInterp();
    const auto& rhoInterp = td.rhoInterp();

    label parceli = 0;
    for (const parcelType& p : c)
    {
        const auto& coords = p.coordinates();
        const auto& tetIs = p.currentTetIndices();

        const vector Uc(UInterp.interpolate(coords, tetIs));

        const scalar pc =
            max
            (
                pInterp.interpolate(coords, tetIs),
                c.constProps().pMin()
            );

        const scalar rhoc(rhoInterp.interpolate(coords, tetIs));
        const scalarField X(liquids.X(p.YLiquid()));
        const scalar sigma = liquids.sigma(pc, p.T(), X);

        We[parceli++] = rhoc*magSqr(p.U() - Uc)*p.d()/sigma;
    }

    const bool haveParticles = c.size();
    if (c.time().writeTime() && returnReduceOr(haveParticles))
    {
        We.write(haveParticles);
    }
}


// ************************************************************************* //
