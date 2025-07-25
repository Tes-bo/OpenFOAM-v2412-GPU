/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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
    Test-ensightFileName

Description
    Check cleanup of ensight file and variable names

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "ensightFileName.H"
#include "ensightVarName.H"
#include "IOstreams.H"

using namespace Foam;

void printCleaning(const fileName& pathName)
{
    Info<< "input = " << pathName << nl;
    Info<< "file  = " << ensight::FileName(pathName) << nl;
    Info<< "var   = " << ensight::VarName(pathName)  << nl;
    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addArgument("fileName .. fileNameN");

    argList args(argc, argv, false, true);

    if (args.size() <= 1 && args.options().empty())
    {
        args.printUsage();
    }

    fileName pathName;

    for (label argI=1; argI < args.size(); ++argI)
    {
        pathName = args[argI];
        printCleaning(pathName);
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
