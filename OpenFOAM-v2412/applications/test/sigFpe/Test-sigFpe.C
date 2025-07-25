/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 Bernhard Gschaider
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
    Test-sigFpe

Description
    Test handling of floating point exceptions by provoking them

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "IOstreams.H"
#include "scalar.H"
#include "sigFpe.H"
#include "argList.H"

#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption("fill-nan", "Test filling memory with NaN");
    argList args(argc, argv);

    // Force on
    sigFpe::unset();
    setEnv("FOAM_SIGFPE", "true", true);
    // setEnv("FOAM_SETNAN", "true", true);

    sigFpe::set(true);

    if (args.found("fill-nan"))
    {
        Info<< nl << "Checking filling with NaN" << endl;
        scalar* data = new scalar[10];

        scalar first = data[0];

        Info<< "First element " << first << endl;
        Info<< "First element times two " << 2*first << endl;

        delete[] data;
    }
    else
    {
        Info<< nl << "Provoking sigFpe (division by zero)" << nl << endl;

        // Writing 1./0. might be optimized away by the compiler

        const scalar zeroVal = ::sin(0.0);
        Info << "Infinity " << 1./zeroVal << endl;
    }

    return 0;
}


// ************************************************************************* //
