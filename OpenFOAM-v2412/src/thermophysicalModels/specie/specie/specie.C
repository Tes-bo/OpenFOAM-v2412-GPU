/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "specie.H"
#include "constants.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(specie, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie::specie(const dictionary& dict)
:
    name_(dict.dictName()),
    Y_(dict.subDict("specie").getOrDefault<scalar>("massFraction", 1)),
    molWeight_(dict.subDict("specie").get<scalar>("molWeight"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::specie::write(Ostream& os) const
{
    // Entries in dictionary format
    {
        os.beginBlock("specie");
        os.writeEntryIfDifferent<scalar>("massFraction", 1, Y_);
        os.writeEntry("molWeight", molWeight_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie& st)
{
    st.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
