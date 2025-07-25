/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "absorptionCoeffs.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::absorptionCoeffs::checkT(const scalar T) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "using absorptionCoeffs out of temperature range:" << nl
            << "    " << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << nl << endl;
    }
}


const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::coeffs
(
    const scalar T
) const
{
    checkT(T);

    if (T < Tcommon_)
    {
        return lowACoeffs_;
    }
    else
    {
        return highACoeffs_;
    }
}


void Foam::radiation::absorptionCoeffs::initialise(const dictionary& dict)
{
    dict.readEntry("Tcommon", Tcommon_);
    dict.readEntry("Tlow", Tlow_);
    dict.readEntry("Thigh", Thigh_);
    dict.readEntry("invTemp", invTemp_);
    dict.readEntry("loTcoeffs", lowACoeffs_);
    dict.readEntry("hiTcoeffs", highACoeffs_);
}


// ************************************************************************* //
