/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "mapDistribute.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistribute::mapDistribute(const dictionary& dict, const label comm)
:
    mapDistribute(comm)
{
    mapDistribute::readDict(dict);
}


Foam::mapDistribute::mapDistribute(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::mapDistribute::readDict(const dictionary& dict)
{
    mapDistributeBase::readDict(dict);

    // Treat missing as empty
    const dictionary* subdict = dict.findDict("transforms");

    if (subdict)
    {
        subdict->readIfPresent("elements", transformElements_);
        subdict->readIfPresent("starts", transformStart_);
    }
    else
    {
        transformElements_.clear();
        transformStart_.clear();
    }
}


void Foam::mapDistribute::writeEntries(Ostream& os) const
{
    mapDistributeBase::writeEntries(os);

    if (transformElements_.size())
    {
        os << nl;
        os.beginBlock("transforms");
        os.writeEntry("elements", transformElements_);
        transformStart_.writeEntry("starts", os);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, mapDistribute& map)
{
    is.fatalCheck(FUNCTION_NAME);

    is  >> static_cast<mapDistributeBase&>(map)
        >> map.transformElements_ >> map.transformStart_;

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const mapDistribute& map)
{
    os  << static_cast<const mapDistributeBase&>(map) << token::NL
        << map.transformElements_ << token::NL
        << map.transformStart_;

    return os;
}


// ************************************************************************* //
