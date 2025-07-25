/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
------------------------------------------------------------------------------
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

#include "faOption.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fa::option::writeHeader(Ostream& os) const
{
    os.beginBlock(name_);
}


void Foam::fa::option::writeFooter(Ostream& os) const
{
    os.endBlock();
}


void Foam::fa::option::writeData(Ostream& os) const
{
    os.writeEntry("type", type());
    os.writeEntry("active", active_);

    os << nl;
    coeffs_.writeEntry(word(type() + "Coeffs"), os);
}


bool Foam::fa::option::read(const dictionary& dict)
{
    dict.readIfPresent("active", active_);

    coeffs_ = dict.optionalSubDict(modelType_ + "Coeffs");

    return true;
}


// ************************************************************************* //
