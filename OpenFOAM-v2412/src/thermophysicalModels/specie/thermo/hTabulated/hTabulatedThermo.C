/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
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

#include "hTabulatedThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hTabulatedThermo<EquationOfState>::hTabulatedThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(dict.subDict("thermodynamics").get<scalar>("Hf")),
    Sf_(dict.subDict("thermodynamics").get<scalar>("Sf")),
    Cp_("Cp", dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hTabulatedThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("thermodynamics");
        os.writeEntry("Hf", Hf_);
        os.writeEntry("Sf", Sf_);
        os.writeEntry("Cp", Cp_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hTabulatedThermo<EquationOfState>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
