/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "ParticleDose.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleDose<CloudType>::ParticleDose
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    GName_(this->coeffDict().getWord("GName"))
{}


template<class CloudType>
Foam::ParticleDose<CloudType>::ParticleDose
(
    const ParticleDose<CloudType>& re
)
:
    CloudFunctionObject<CloudType>(re),
    GName_(re.GName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleDose<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    auto& c = this->owner();

    auto* resultPtr = c.template getObjectPtr<IOField<scalar>>("D");

    if (!resultPtr)
    {
        resultPtr = new IOField<scalar>
        (
            IOobject
            (
                "D",
                c.time().timeName(),
                c,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        );

        resultPtr->store();
    }
    auto& D = *resultPtr;

    D.resize(c.size(), Zero);

    const fvMesh& mesh = this->owner().mesh();

    const auto& G = mesh.lookupObject<volScalarField>(GName_);

    label parceli = 0;
    for (const parcelType& p : c)
    {
        D[parceli] += G[p.cell()]*mesh.time().deltaTValue();
        parceli++;
    }

    const bool haveParticles = c.size();
    if (c.time().writeTime() && returnReduceOr(haveParticles))
    {
        D.write(haveParticles);
    }
}


// ************************************************************************* //
