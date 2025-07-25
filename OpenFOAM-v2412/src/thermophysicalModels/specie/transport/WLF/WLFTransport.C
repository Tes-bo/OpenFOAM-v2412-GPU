/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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

#include "WLFTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::WLFTransport<Thermo>::WLFTransport(const dictionary& dict)
:
    Thermo(dict),
    mu0_(readCoeff("mu0", dict)),
    Tr_(readCoeff("Tr", dict)),
    C1_(readCoeff("C1", dict)),
    C2_(readCoeff("C2", dict)),
    rPr_(1.0/readCoeff("Pr", dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::WLFTransport<Thermo>::write(Ostream& os) const
{
    os.beginBlock(this->specie::name());

    Thermo::write(os);

    {
        os.beginBlock("transport");
        os.writeEntry("mu0", mu0_);
        os.writeEntry("Tr", Tr_);
        os.writeEntry("C1", C1_);
        os.writeEntry("C2", C2_);
        os.writeEntry("Pr", 1.0/rPr_);
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const WLFTransport<Thermo>& wlft
)
{
    wlft.write(os);
    return os;
}


// ************************************************************************* //
