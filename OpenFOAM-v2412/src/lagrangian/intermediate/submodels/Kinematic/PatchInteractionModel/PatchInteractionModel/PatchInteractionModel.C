/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "PatchInteractionModel.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::PatchInteractionModel<CloudType>::interactionTypeNames_
{
    "rebound", "stick", "escape"
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchInteractionModel<CloudType>::writeFileHeader(Ostream& os)
{
    this->writeHeader(os, "Particle patch interaction");
    this->writeHeaderValue(os, "Model", this->modelType());

    this->writeCommented(os, "Time");
    this->writeTabbed(os, "escapedParcels");
    this->writeTabbed(os, "escapedMass");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::word Foam::PatchInteractionModel<CloudType>::interactionTypeToWord
(
    const interactionType& itEnum
)
{
    word it = "other";

    switch (itEnum)
    {
        case itNone:
        {
            it = "none";
            break;
        }
        case itRebound:
        {
            it = "rebound";
            break;
        }
        case itStick:
        {
            it = "stick";
            break;
        }
        case itEscape:
        {
            it = "escape";
            break;
        }
        default:
        {
        }
    }

    return it;
}


template<class CloudType>
typename Foam::PatchInteractionModel<CloudType>::interactionType
Foam::PatchInteractionModel<CloudType>::wordToInteractionType
(
    const word& itWord
)
{
    if (itWord == "none")
    {
        return itNone;
    }
    if (itWord == "rebound")
    {
        return itRebound;
    }
    else if (itWord == "stick")
    {
        return itStick;
    }
    else if (itWord == "escape")
    {
        return itEscape;
    }
    else
    {
        return itOther;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner),
    functionObjects::writeFile(owner, this->localPath(), typeName, false),
    UName_("unknown_U"),
    escapedParcels_(0),
    escapedMass_(0.0),
    Urmax_(1e-4)
{}


template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    functionObjects::writeFile
    (
        owner,
        this->localPath(),
        type,
        this->coeffDict(),
        false                   // Do not write by default
    ),
    UName_(this->coeffDict().template getOrDefault<word>("U", "U")),
    escapedParcels_(0),
    escapedMass_(0.0),
    Urmax_(this->coeffDict().template getOrDefault<scalar>("UrMax", 0))
{}


template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    const PatchInteractionModel<CloudType>& pim
)
:
    CloudSubModelBase<CloudType>(pim),
    functionObjects::writeFile(pim),
    UName_(pim.UName_),
    escapedParcels_(pim.escapedParcels_),
    escapedMass_(pim.escapedMass_),
    Urmax_(pim.Urmax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::word& Foam::PatchInteractionModel<CloudType>::UName() const
{
    return UName_;
}


template<class CloudType>
const Foam::scalar& Foam::PatchInteractionModel<CloudType>::Urmax() const
{
    return Urmax_;
}


template<class CloudType>
void Foam::PatchInteractionModel<CloudType>::addToEscapedParcels
(
    const scalar mass
)
{
    escapedMass_ += mass;
    ++escapedParcels_;
}


template<class CloudType>
void Foam::PatchInteractionModel<CloudType>::postEvolve()
{}


template<class CloudType>
void Foam::PatchInteractionModel<CloudType>::info()
{
    CloudSubModelBase<CloudType>::info();

    const label escapedParcels0 =
        this->template getBaseProperty<label>("escapedParcels");
    const label escapedParcelsTotal =
        escapedParcels0 + returnReduce(escapedParcels_, sumOp<label>());

    const scalar escapedMass0 =
        this->template getBaseProperty<scalar>("escapedMass");
    const scalar escapedMassTotal =
        escapedMass0 + returnReduce(escapedMass_, sumOp<scalar>());

    Log_<< "    Parcel fate: system (number, mass)" << nl
        << "      - escape                      = " << escapedParcelsTotal
        << ", " << escapedMassTotal << endl;

    if (!this->writtenHeader_)
    {
        this->writeFileHeader(this->file());
        this->writtenHeader_ = true;
        this->file() << endl;
    }

    this->writeCurrentTime(this->file());
    this->file()
        << tab << escapedParcelsTotal << tab << escapedMassTotal;


    if (this->writeTime())
    {
        this->setBaseProperty("escapedParcels", escapedParcelsTotal);
        escapedParcels_ = 0;

        this->setBaseProperty("escapedMass", escapedMassTotal);
        escapedMass_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchInteractionModelNew.C"

// ************************************************************************* //
