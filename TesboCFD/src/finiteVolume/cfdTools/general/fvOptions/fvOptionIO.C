/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "fvOption.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::option::writeHeader(Ostream& os) const
{
    os.beginBlock(name_);
}


void Foam::fv::option::writeFooter(Ostream& os) const
{
    os.endBlock();
}


void Foam::fv::option::writeData(Ostream& os) const
{
    os.writeEntry("type", type());
    os.writeEntry("active", active_);

    os << nl;
    coeffs_.writeEntry(word(type() + "Coeffs"), os);
}


bool Foam::fv::option::read(const dictionary& dict)
{
    dict.readIfPresent("active", active_);

    log = dict.getOrDefault("log", true);

    coeffs_ = dict.optionalSubDict(modelType_ + "Coeffs");

    return true;
}


// ************************************************************************* //
