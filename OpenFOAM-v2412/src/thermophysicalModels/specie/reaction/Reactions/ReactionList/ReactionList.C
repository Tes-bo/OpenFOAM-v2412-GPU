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

#include "ReactionList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const ReactionTable<ThermoType>& thermoDb
)
:
    SLPtrList<Reaction<ThermoType>>(),
    species_(species),
    thermoDb_(thermoDb),
    dict_()
{}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const ReactionTable<ThermoType>& thermoDb,
    const dictionary& dict
)
:
    SLPtrList<Reaction<ThermoType>>(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dict)
{
    readReactionDict();
}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList(const ReactionList& reactions)
:
    SLPtrList<Reaction<ThermoType>>(reactions),
    species_(reactions.species_),
    thermoDb_(reactions.thermoDb_),
    dict_(reactions.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::ReactionList<ThermoType>::readReactionDict()
{
    for (const entry& dEntry : dict_.subDict("reactions"))
    {
        this->append
        (
            Reaction<ThermoType>::New
            (
                species_,
                thermoDb_,
                dEntry.dict()
            ).ptr()
        );
    }

    return true;
}


template<class ThermoType>
void Foam::ReactionList<ThermoType>::write(Ostream& os) const
{
    os.beginBlock("reactions");

    for (const Reaction<ThermoType>& r : *this)
    {
        os.beginBlock(r.name());

        os.writeEntry("type", r.type());
        r.write(os);

        os.endBlock();
    }

    os.endBlock();
}


// ************************************************************************* //
