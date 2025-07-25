/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "CloudFunctionObject.H"
#include "functionObject.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::write()
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    outputDir_()
{}


template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName,
    const word& objectType
)
:
    CloudSubModelBase<CloudType>(modelName, owner, dict, typeName, objectType),
    outputDir_()
{
    // Put in undecomposed case
    // (Note: gives problems for distributed data running)

    outputDir_ =
    (
        owner.mesh().time().globalPath()
      / functionObject::outputPrefix
      / this->localPath()
    );

    outputDir_.clean();  // Remove unneeded ".."
}


template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject
(
    const CloudFunctionObject<CloudType>& ppm
)
:
    CloudSubModelBase<CloudType>(ppm),
    outputDir_(ppm.outputDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    if (this->owner().time().writeTime())
    {
        this->write();
    }
}


template<class CloudType>
bool Foam::CloudFunctionObject<CloudType>::postMove
(
    parcelType&,
    const scalar,
    const point&,
    const typename parcelType::trackingData& td
)
{
    return true;
}


template<class CloudType>
bool Foam::CloudFunctionObject<CloudType>::postPatch
(
    const parcelType&,
    const polyPatch&,
    const typename parcelType::trackingData& td
)
{
    return true;
}


template<class CloudType>
bool Foam::CloudFunctionObject<CloudType>::postFace
(
    const parcelType&,
    const typename parcelType::trackingData& td
)
{
    return true;
}


template<class CloudType>
const Foam::fileName& Foam::CloudFunctionObject<CloudType>::outputDir() const
{
    return outputDir_;
}


template<class CloudType>
Foam::fileName Foam::CloudFunctionObject<CloudType>::writeTimeDir() const
{
    return outputDir_/this->owner().time().timeName();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CloudFunctionObjectNew.C"

// ************************************************************************* //
