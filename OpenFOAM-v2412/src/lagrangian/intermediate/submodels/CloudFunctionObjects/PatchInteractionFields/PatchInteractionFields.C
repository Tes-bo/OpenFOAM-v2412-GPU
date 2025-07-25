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

#include "PatchInteractionFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::Enum<typename Foam::PatchInteractionFields<CloudType>::resetMode>
Foam::PatchInteractionFields<CloudType>::resetModeNames_
({
    { resetMode::none, "none" },
    { resetMode::timeStep, "timeStep" },
    { resetMode::writeTime, "writeTime" },
});

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchInteractionFields<CloudType>::clearOrReset
(
    autoPtr<volScalarField>& fieldPtr,
    const word& fieldName,
    const dimensionSet& dims
) const
{
    if (fieldPtr)
    {
        fieldPtr->primitiveFieldRef() = scalar(0);
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        fieldPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::scopedName
                    (
                        this->owner().name(),
                        this->modelName(),
                        fieldName
                    ),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dims, Zero)
            )
        );
    }
}


template<class CloudType>
void Foam::PatchInteractionFields<CloudType>::reset()
{
    clearOrReset(massPtr_, "mass", dimMass);
    clearOrReset(countPtr_, "count", dimless);
}


template<class CloudType>
void Foam::PatchInteractionFields<CloudType>::write()
{
    if (massPtr_)
    {
        massPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "massPtr not valid" << abort(FatalError);
    }

    if (countPtr_)
    {
        countPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "countPtr not valid" << abort(FatalError);
    }

    if (resetMode_ == resetMode::writeTime)
    {
        reset();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInteractionFields<CloudType>::PatchInteractionFields
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    massPtr_(nullptr),
    countPtr_(nullptr),
    resetMode_
    (
        resetModeNames_.getOrDefault
        (
            "resetMode",
            this->coeffDict(),
            resetMode::none
        )
    )
{
    reset();
}


template<class CloudType>
Foam::PatchInteractionFields<CloudType>::PatchInteractionFields
(
    const PatchInteractionFields<CloudType>& pii
)
:
    CloudFunctionObject<CloudType>(pii),
    massPtr_(nullptr),
    countPtr_(nullptr),
    resetMode_(pii.resetMode_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchInteractionFields<CloudType>::preEvolve
(
    const typename parcelType::trackingData&
)
{
    if (resetMode_ == resetMode::timeStep)
    {
        reset();
    }
}


template<class CloudType>
bool Foam::PatchInteractionFields<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    const typename parcelType::trackingData& td
)
{
    const label patchi = pp.index();
    const label facei = pp.whichFace(p.face());

    massPtr_->boundaryFieldRef()[patchi][facei] += p.nParticle()*p.mass();
    countPtr_->boundaryFieldRef()[patchi][facei] += 1;

    return true;
}


// ************************************************************************* //
