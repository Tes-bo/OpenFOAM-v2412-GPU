/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Test-io

Description
    Test basic stream functionality

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOmanip.H"
#include "scalar.H"
#include "List.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(void)
{
    string st("sfdsf  sdfs23df sdf32f .  sdfsdff23/2sf32");
    Info<<"string: " << st << nl;
    Info<<"word:   \"" << string::validate<word>(st) << "\"" << endl;

    string st1("1234567");

    Info<< label(st1.size()) << tab << string(word(st1)) << endl;

    Info<< setw(20) << setprecision(3) << 1.234234 << endl;

    Info<< hex << 255 << endl;


    Info<< nl << "Formatted fields" << nl;

    const char oldfill = static_cast<OSstream&>(Info).fill();

    Info<< setfill('-');
    Info<< "|" << setf(ios_base::left) << setw(32) << " foo " << "|" << nl;
    Info<< "|" << setf(ios_base::left) << setw(10) << " bar " << "|" << nl;
    Info<< "|" << setf(ios_base::left) << setw(10) << "" << "|" << nl;

    Info<< "resetting fill from (0x" << hex << int(oldfill) << ")" << nl;
    Info<< setfill(oldfill);
    Info<< "|" << setf(ios_base::left) << setw(10) << " bar " << "|" << nl;

    Info<< nl << nl;

    Info.operator Foam::OSstream&() << "stop" << endl;

    static_cast<OSstream&>(Info) << "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
