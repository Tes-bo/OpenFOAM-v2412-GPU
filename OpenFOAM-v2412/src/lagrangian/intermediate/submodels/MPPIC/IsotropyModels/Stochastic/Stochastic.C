/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "Stochastic.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropyModels::Stochastic<CloudType>::Stochastic
(
    const dictionary& dict,
    CloudType& owner
)
:
    IsotropyModel<CloudType>(dict, owner, typeName)
{}


template<class CloudType>
Foam::IsotropyModels::Stochastic<CloudType>::Stochastic
(
    const Stochastic<CloudType>& cm
)
:
    IsotropyModel<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropyModels::Stochastic<CloudType>::~Stochastic()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::IsotropyModels::Stochastic<CloudType>::sampleGauss()
{
    static bool isCached = true;
    static scalar xCached;

    if (isCached)
    {
        isCached = false;

        return xCached;
    }
    else
    {
        Random& rndGen = this->owner().rndGen();

        scalar f, m, x, y;

        do
        {
            x = 2.0*rndGen.sample01<scalar>() - 1.0;
            y = 2.0*rndGen.sample01<scalar>() - 1.0;
            m = x*x + y*y;
        } while (m >= 1.0 || m == 0.0);

        f = sqrt(-2.0*log(m)/m);
        xCached = x*f;
        isCached = true;

        return y*f;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::IsotropyModels::Stochastic<CloudType>::calculate()
{
    const fvMesh& mesh = this->owner().mesh();
    const scalar deltaT(this->owner().db().time().deltaTValue());
    Random& rndGen = this->owner().rndGen();

    const scalar oneBySqrtThree = sqrt(1.0/3.0);

    const auto& volumeAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            IOobject::scopedName(this->owner().name(), "volumeAverage")
        );
    const auto& radiusAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            IOobject::scopedName(this->owner().name(), "radiusAverage")
        );
    const auto& uAverage =
        mesh.lookupObject<AveragingMethod<vector>>
        (
            IOobject::scopedName(this->owner().name(), "uAverage")
        );
    const auto& uSqrAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            IOobject::scopedName(this->owner().name(), "uSqrAverage")
        );
    const auto& frequencyAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            IOobject::scopedName(this->owner().name(), "frequencyAverage")
        );
    const auto& massAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            IOobject::scopedName(this->owner().name(), "massAverage")
        );

    // calculate time scales and pdf exponent
    autoPtr<AveragingMethod<scalar>> exponentAveragePtr
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(this->owner().name(), "exponentAverage"),
                this->owner().db().time().timeName(),
                mesh
            ),
            this->owner().solution().dict(),
            mesh
        )
    );
    AveragingMethod<scalar>& exponentAverage = exponentAveragePtr();
    exponentAverage =
        exp
        (
          - deltaT
           *this->timeScaleModel_->oneByTau
            (
                volumeAverage,
                radiusAverage,
                uSqrAverage,
                frequencyAverage
            )
        )();

    // random sampling
    for (typename CloudType::parcelType& p : this->owner())
    {
        const tetIndices tetIs(p.currentTetIndices());

        const scalar x = exponentAverage.interpolate(p.coordinates(), tetIs);

        if (x < rndGen.sample01<scalar>())
        {
            const vector r(sampleGauss(), sampleGauss(), sampleGauss());

            const vector u = uAverage.interpolate(p.coordinates(), tetIs);
            const scalar uRms =
                sqrt(max(uSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));

            p.U() = u + r*uRms*oneBySqrtThree;
        }
    }

    // correction velocity averages
    autoPtr<AveragingMethod<vector>> uTildeAveragePtr
    (
        AveragingMethod<vector>::New
        (
            IOobject
            (
                IOobject::scopedName(this->owner().name(), "uTildeAverage"),
                this->owner().db().time().timeName(),
                mesh
            ),
            this->owner().solution().dict(),
            mesh
        )
    );
    AveragingMethod<vector>& uTildeAverage = uTildeAveragePtr();
    for (typename CloudType::parcelType& p : this->owner())
    {
        const tetIndices tetIs(p.currentTetIndices());
        uTildeAverage.add(p.coordinates(), tetIs, p.nParticle()*p.mass()*p.U());
    }
    uTildeAverage.average(massAverage);

    autoPtr<AveragingMethod<scalar>> uTildeSqrAveragePtr
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(this->owner().name(), "uTildeSqrAverage"),
                this->owner().db().time().timeName(),
                mesh
            ),
            this->owner().solution().dict(),
            mesh
        )
    );
    AveragingMethod<scalar>& uTildeSqrAverage = uTildeSqrAveragePtr();
    for (typename CloudType::parcelType& p : this->owner())
    {
        const tetIndices tetIs(p.currentTetIndices());
        const vector uTilde = uTildeAverage.interpolate(p.coordinates(), tetIs);
        uTildeSqrAverage.add
        (
            p.coordinates(),
            tetIs,
            p.nParticle()*p.mass()*magSqr(p.U() - uTilde)
        );
    }
    uTildeSqrAverage.average(massAverage);

    // conservation correction
    for (typename CloudType::parcelType& p : this->owner())
    {
        const tetIndices tetIs(p.currentTetIndices());

        const vector u = uAverage.interpolate(p.coordinates(), tetIs);
        const scalar uRms =
            sqrt(max(uSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));

        const vector uTilde = uTildeAverage.interpolate(p.coordinates(), tetIs);
        const scalar uTildeRms =
            sqrt
            (
                max(uTildeSqrAverage.interpolate(p.coordinates(), tetIs), 0.0)
            );

        p.U() = u + (p.U() - uTilde)*uRms/max(uTildeRms, SMALL);
    }
}


// ************************************************************************* //
