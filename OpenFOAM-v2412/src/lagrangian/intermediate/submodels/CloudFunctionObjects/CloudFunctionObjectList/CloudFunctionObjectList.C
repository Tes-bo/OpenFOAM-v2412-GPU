/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "CloudFunctionObjectList.H"
#include "entry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObjectList<CloudType>::CloudFunctionObjectList
(
    CloudType& owner
)
:
    PtrList<CloudFunctionObject<CloudType>>(),
    owner_(owner),
    dict_()
{}


template<class CloudType>
Foam::CloudFunctionObjectList<CloudType>::CloudFunctionObjectList
(
    CloudType& owner,
    const dictionary& dict,
    const bool readFields
)
:
    PtrList<CloudFunctionObject<CloudType>>(),
    owner_(owner),
    dict_(dict)
{
    if (readFields)
    {
        Info<< "Constructing cloud functions" << endl;

        this->resize(dict.size());

        label count = 0;
        for (const word& modelName : dict.toc())
        {
            const dictionary& modelDict = dict.subDict(modelName);

            {
                this->set
                (
                    count,
                    CloudFunctionObject<CloudType>::New
                    (
                        modelDict,
                        owner,
                        modelDict.get<word>("type"),
                        modelName
                    )
                );
            }
            ++count;
        }

        if (!count)
        {
            Info<< "    none" << endl;
        }
    }
}


template<class CloudType>
Foam::CloudFunctionObjectList<CloudType>::CloudFunctionObjectList
(
    const CloudFunctionObjectList& cfol
)
:
    PtrList<CloudFunctionObject<CloudType>>(cfol),
    owner_(cfol.owner_),
    dict_(cfol.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObjectList<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    for (auto& cfo : *this)
    {
        cfo.preEvolve(td);
    }
}


template<class CloudType>
void Foam::CloudFunctionObjectList<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    for (auto& cfo : *this)
    {
        cfo.postEvolve(td);
    }
}


template<class CloudType>
bool Foam::CloudFunctionObjectList<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point& position0,
    const typename parcelType::trackingData& td
)
{
    if (!td.keepParticle)
    {
        return td.keepParticle;
    }

    for (auto& cfo : *this)
    {
        if (!cfo.postMove(p, dt, position0, td))
        {
            return false;
        }
    }

    return true;
}


template<class CloudType>
bool Foam::CloudFunctionObjectList<CloudType>::postPatch
(
    parcelType& p,
    const polyPatch& pp,
    const typename parcelType::trackingData& td
)
{
    if (!td.keepParticle)
    {
        return td.keepParticle;
    }

    for (auto& cfo : *this)
    {
        if (!cfo.postPatch(p, pp, td))
        {
            return false;
        }
    }

    return true;
}


template<class CloudType>
bool Foam::CloudFunctionObjectList<CloudType>::postFace
(
    parcelType& p,
    const typename parcelType::trackingData& td
)
{
    if (!td.keepParticle)
    {
        return td.keepParticle;
    }

    for (auto& cfo : *this)
    {
        if (!cfo.postFace(p, td))
        {
            return false;
        }
    }

    return true;
}


// ************************************************************************* //
