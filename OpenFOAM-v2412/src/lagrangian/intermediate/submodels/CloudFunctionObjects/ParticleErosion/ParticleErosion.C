/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "ParticleErosion.H"
#include "wordRes.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleErosion<CloudType>::resetQ()
{
    if (QPtr_)
    {
        QPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        QPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Q",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimVolume, Zero)
            )
        );
    }
}


template<class CloudType>
Foam::label Foam::ParticleErosion<CloudType>::applyToPatch
(
    const label globalPatchi
) const
{
    return patchIDs_.find(globalPatchi);
}


template<class CloudType>
void Foam::ParticleErosion<CloudType>::write()
{
    if (QPtr_)
    {
        QPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "QPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleErosion<CloudType>::ParticleErosion
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    QPtr_(nullptr),
    patchIDs_(),
    p_(this->coeffDict().getScalar("p")),
    psi_(this->coeffDict().template getOrDefault<scalar>("psi", 2.0)),
    K_(this->coeffDict().template getOrDefault<scalar>("K", 2.0))
{
    const wordList allPatchNames(owner.mesh().boundaryMesh().names());
    const wordRes patchNames
    (
        this->coeffDict().template get<wordRes>("patches")
    );

    labelHashSet uniqIds;
    for (const wordRe& select : patchNames)
    {
        labelList ids = wordRes::matching(select, allPatchNames);

        if (ids.empty())
        {
            WarningInFunction
                << "Cannot find any patch names matching "
                << select << nl;
        }

        uniqIds.insert(ids);
    }

    patchIDs_ = uniqIds.sortedToc();

    // Trigger creation of the Q field
    resetQ();
}


template<class CloudType>
Foam::ParticleErosion<CloudType>::ParticleErosion
(
    const ParticleErosion<CloudType>& pe
)
:
    CloudFunctionObject<CloudType>(pe),
    QPtr_(nullptr),
    patchIDs_(pe.patchIDs_),
    p_(pe.p_),
    psi_(pe.psi_),
    K_(pe.K_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleErosion<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    resetQ();
}


template<class CloudType>
bool Foam::ParticleErosion<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    const typename parcelType::trackingData& td
)
{
    const label patchi = pp.index();

    const label localPatchi = applyToPatch(patchi);

    if (localPatchi != -1)
    {
        vector nw;
        vector Up;

        // patch-normal direction
        this->owner().patchData(p, pp, nw, Up);

        // particle velocity relative to patch
        const vector& U = p.U() - Up;

        // quick reject if particle travelling away from the patch
        if ((nw & U) < 0)
        {
            return true;
        }

        const scalar magU = mag(U);
        const vector Udir = U/magU;

        // determine impact angle, alpha
        const scalar alpha = mathematical::piByTwo - acos(nw & Udir);

        const scalar coeff = p.nParticle()*p.mass()*sqr(magU)/(p_*psi_*K_);

        const label patchFacei = pp.whichFace(p.face());
        scalar& Q = QPtr_->boundaryFieldRef()[patchi][patchFacei];
        if (tan(alpha) < K_/6.0)
        {
            Q += coeff*(sin(2.0*alpha) - 6.0/K_*sqr(sin(alpha)));
        }
        else
        {
            Q += coeff*(K_*sqr(cos(alpha))/6.0);
        }
    }

    return true;
}


// ************************************************************************* //
