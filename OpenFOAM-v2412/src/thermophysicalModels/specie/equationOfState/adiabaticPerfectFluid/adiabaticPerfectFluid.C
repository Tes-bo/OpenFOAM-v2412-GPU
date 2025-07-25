/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "adiabaticPerfectFluid.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const dictionary& dict
)
:
    Specie(dict),
    p0_(dict.subDict("equationOfState").get<scalar>("p0")),
    rho0_(dict.subDict("equationOfState").get<scalar>("rho0")),
    gamma_(dict.subDict("equationOfState").get<scalar>("gamma")),
    B_(dict.subDict("equationOfState").get<scalar>("B"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::adiabaticPerfectFluid<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("equationOfState");
        os.writeEntry("p0", p0_);
        os.writeEntry("rho0", rho0_);
        os.writeEntry("gamma", gamma_);
        os.writeEntry("B", B_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const adiabaticPerfectFluid<Specie>& pf
)
{
    pf.write(os);
    return os;
}


// ************************************************************************* //
