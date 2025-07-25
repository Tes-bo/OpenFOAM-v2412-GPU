/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

Application
    Test-unitConversion

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "unitConversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<< "30_deg:   " << 30_deg << nl;
    Info<< "30.0_deg: " << 30.0_deg << nl;
    Info<< "3e+1_deg: " << 3e+1_deg << nl;
    Info<< "degToRad(30): " << degToRad(30) << nl;

    Info<< "cos(30_deg): " << ::cos(30_deg) << nl;
    Info<< "1000 rpm = " << rpmToRads(1000) << " 1/s" << nl;
    Info<< "100 1/s  = " << radsToRpm(100) << " rpm" << nl;

    return 0;
}


// ************************************************************************* //
