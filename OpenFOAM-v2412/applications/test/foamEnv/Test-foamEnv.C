/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application
    Test-etcFiles

Description
    Test etcFiles functionality.
    Similar to foamEtcFile script, but automatically prunes nonexistent
    directories from the list.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::removeOption("case");

    argList::addArgument("env...");

    argList::addNote
    (
        "Simple test/report OpenFOAM environment"
    );

    argList args(argc, argv, false, true);

    for (int argi = 1; argi < args.size(); ++argi)
    {
        const std::string envName(args[argi]);

        if (hasEnv(envName))
        {
            Info<<"Have env " << envName.c_str() << "=" << getEnv(envName)
                << nl;
        }
        else
        {
            Info<<"No env " << envName.c_str()<< nl;
        }
    }

    return 0;
}


// ************************************************************************* //
