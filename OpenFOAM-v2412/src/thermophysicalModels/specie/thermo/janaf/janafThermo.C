/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "janafThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janafThermo<EquationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorInFunction
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::janafThermo<EquationOfState>::janafThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Tlow_(dict.subDict("thermodynamics").get<scalar>("Tlow")),
    Thigh_(dict.subDict("thermodynamics").get<scalar>("Thigh")),
    Tcommon_(dict.subDict("thermodynamics").get<scalar>("Tcommon")),
    highCpCoeffs_(dict.subDict("thermodynamics").lookup("highCpCoeffs")),
    lowCpCoeffs_(dict.subDict("thermodynamics").lookup("lowCpCoeffs"))
{
    // Convert coefficients to mass-basis
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] *= this->R();
        lowCpCoeffs_[coefLabel] *= this->R();
    }

    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janafThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    // Convert coefficients back to dimensionless form
    coeffArray highCpCoeffs;
    coeffArray lowCpCoeffs;
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs[coefLabel] = highCpCoeffs_[coefLabel]/this->R();
        lowCpCoeffs[coefLabel] = lowCpCoeffs_[coefLabel]/this->R();
    }

    // Entries in dictionary format
    {
        os.beginBlock("thermodynamics");
        os.writeEntry("Tlow", Tlow_);
        os.writeEntry("Thigh", Thigh_);
        os.writeEntry("Tcommon", Tcommon_);
        os.writeEntry("highCpCoeffs", highCpCoeffs);
        os.writeEntry("lowCpCoeffs", lowCpCoeffs);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const janafThermo<EquationOfState>& jt
)
{
    jt.write(os);
    return os;
}


// ************************************************************************* //
