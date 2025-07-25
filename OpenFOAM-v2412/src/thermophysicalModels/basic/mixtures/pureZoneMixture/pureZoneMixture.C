/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021,2022 OpenCFD Ltd.
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

#include "pureZoneMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::pureZoneMixture<ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    const auto& czs = mesh_.cellZones();

    const auto* dictPtr = thermoDict.findDict("none");

    speciesData_.setSize(dictPtr ? czs.size()+1 : czs.size());
    forAll(czs, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.subDict(czs[i].name()))
        );
    }

    if (dictPtr)
    {
        speciesData_.set(czs.size(), new ThermoType(*dictPtr));
    }

    return speciesData_[0];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class ThermoType>
Foam::pureZoneMixture<ThermoType>::pureZoneMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMixture(thermoDict, mesh, phaseName),
    mesh_(mesh),
    mixture_("mixture", constructSpeciesData(thermoDict.subDict("mixture")))
{
    // Cache index per cell. This is the cellZone except for unzoned cells
    // which are last

    const auto& czs = mesh_.cellZones();
    zoneID_.setSize(mesh_.nCells(), czs.size());
    for (const auto& cz : czs)
    {
        UIndirectList<label>(zoneID_, cz) = cz.index();
    }

    // Check if any unzoned cells but no 'none' specification
    if (speciesData_.size() == czs.size())
    {
        const label noneCelli = zoneID_.find(czs.size());
        if (noneCelli != -1)
        {
            FatalErrorInFunction << "Have unzoned cell " << noneCelli
                << " at " << mesh_.cellCentres()[noneCelli]
                << " but no \"none\" entry in \"mixture\""
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::pureZoneMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = speciesData_[zoneID_[celli]];
    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::pureZoneMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    const label celli = mesh_.boundary()[patchi].faceCells()[facei];
    mixture_ = speciesData_[zoneID_[celli]];
    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::pureZoneMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    // (per zone) constant density
    return this->cellMixture(celli);
}


template<class ThermoType>
const ThermoType& Foam::pureZoneMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    return this->patchFaceMixture(patchi, facei);
}


template<class ThermoType>
void Foam::pureZoneMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    constructSpeciesData(thermoDict);
}


// ************************************************************************* //
