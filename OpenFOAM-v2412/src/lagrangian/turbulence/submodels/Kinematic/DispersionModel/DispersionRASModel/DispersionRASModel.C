/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "DispersionRASModel.H"
#include "demandDrivenData.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::DispersionRASModel<CloudType>::kModel() const
{
    const objectRegistry& obr = this->owner().mesh();
    const word turbName =
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            this->owner().U().group()
        );

    const turbulenceModel* turb = obr.findObject<turbulenceModel>(turbName);

    if (turb)
    {
        return turb->k();
    }

    FatalErrorInFunction
        << "Turbulence model not found in mesh database" << nl
        << "Database objects include: " << obr.sortedToc()
        << abort(FatalError);

    return nullptr;
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::DispersionRASModel<CloudType>::epsilonModel() const
{
    const objectRegistry& obr = this->owner().mesh();
    const word turbName =
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            this->owner().U().group()
        );

    const turbulenceModel* turb = obr.findObject<turbulenceModel>(turbName);

    if (turb)
    {
        return turb->epsilon();
    }

    FatalErrorInFunction
        << "Turbulence model not found in mesh database" << nl
        << "Database objects include: " << obr.sortedToc()
        << abort(FatalError);

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DispersionRASModel<CloudType>::DispersionRASModel
(
    const dictionary&,
    CloudType& owner
)
:
    DispersionModel<CloudType>(owner),
    kPtr_(nullptr),
    epsilonPtr_(nullptr),
    ownK_(false),
    ownEpsilon_(false)
{}


template<class CloudType>
Foam::DispersionRASModel<CloudType>::DispersionRASModel
(
    const DispersionRASModel<CloudType>& dm
)
:
    DispersionModel<CloudType>(dm),
    kPtr_(dm.kPtr_),
    epsilonPtr_(dm.epsilonPtr_),
    ownK_(dm.ownK_),
    ownEpsilon_(dm.ownEpsilon_)
{
    dm.ownK_ = false;
    dm.ownEpsilon_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DispersionRASModel<CloudType>::~DispersionRASModel()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DispersionRASModel<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        tmp<volScalarField> tk = this->kModel();
        if (tk.movable())
        {
            // Take ownership
            kPtr_ = tk.ptr();
            ownK_ = true;
        }
        else
        {
            kPtr_ = &tk.cref();
            ownK_ = false;
        }

        tmp<volScalarField> tepsilon = this->epsilonModel();
        if (tepsilon.movable())
        {
            // Take ownership
            epsilonPtr_ = tepsilon.ptr();
            ownEpsilon_ = true;
        }
        else
        {
            epsilonPtr_ = &tepsilon();
            ownEpsilon_ = false;
        }
    }
    else
    {
        if (ownK_)
        {
            deleteDemandDrivenData(kPtr_);
            ownK_ = false;
        }
        if (ownEpsilon_)
        {
            deleteDemandDrivenData(epsilonPtr_);
            ownEpsilon_ = false;
        }
    }
}


template<class CloudType>
void Foam::DispersionRASModel<CloudType>::write(Ostream& os) const
{
    DispersionModel<CloudType>::write(os);

    os.writeEntry("ownK", ownK_);
    os.writeEntry("ownEpsilon", ownEpsilon_);
}


// ************************************************************************* //
