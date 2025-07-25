/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addArgument("file .. fileN");

    argList::removeOption("case");
    argList::addOption("ext", "bak");

    argList args(argc, argv, false, true);

    if (args.size() <= 1)
    {
        args.printUsage();
    }

    label ok = 0;

    for (label argi=1; argi < args.size(); ++argi)
    {
        const auto srcFile = args.get<fileName>(argi);

        if (args.found("ext"))
        {
            if (mvBak(srcFile, args["ext"]))
            {
                ok++;
            }
        }
        else
        {
            if (mvBak(srcFile))
            {
                ok++;
            }
        }
    }

    Info<< "mvBak called for " << args.size()-1
        << " files (moved " << ok << ")\n" << endl;

    return 0;
}


// ************************************************************************* //
