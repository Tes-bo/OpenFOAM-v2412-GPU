/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "hPolynomialThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(dict.subDict("thermodynamics").get<scalar>("Hf")),
    Sf_(dict.subDict("thermodynamics").get<scalar>("Sf")),
    CpCoeffs_(dict.subDict("thermodynamics").lookup(coeffsName("Cp"))),
    Tref_(dict.subDict("thermodynamics").getOrDefault<scalar>("Tref", Tstd)),
    Href_(dict.subDict("thermodynamics").getOrDefault<scalar>("Href", 0)),
    Sref_(dict.subDict("thermodynamics").getOrDefault<scalar>("Sref", 0)),
    Pref_(dict.subDict("thermodynamics").getOrDefault<scalar>("Pref", Pstd)),
    hCoeffs_(),
    sCoeffs_()
{
    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();
    // Offset h poly so that it is relative to the enthalpy at Tref
    hCoeffs_[0] +=
        Href_ - hCoeffs_.value(Tref_) - EquationOfState::H(Pstd, Tref_);

    // Offset s poly so that it is relative to the entropy at Tref
    sCoeffs_[0] +=
        Sref_ - sCoeffs_.value(Tref_) - EquationOfState::S(Pstd, Tref_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
void Foam::hPolynomialThermo<EquationOfState, PolySize>::write
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
        os.writeEntry("Tref", Tref_);
        os.writeEntry("Hsref", Href_);
        os.writeEntry("Sref", Sref_);
        os.writeEntry("Pref", Pref_);
        os.writeEntry(coeffsName("Cp"), CpCoeffs_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
