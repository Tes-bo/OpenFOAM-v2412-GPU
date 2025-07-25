/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Test-fieldTypes

Description
    Print fieldTypes

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"

#include "areaFields.H"
#include "fieldTypes.H"
#include "pointFields.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

    Info<< "basic:    " << flatOutput(fieldTypes::basic) << nl
        << "area:     " << flatOutput(fieldTypes::area) << nl
        << "volume:   " << flatOutput(fieldTypes::volume) << nl
        << "internal: " << flatOutput(fieldTypes::internal) << nl
        << "point:    " << flatOutput(fieldTypes::point) << nl
        << endl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
