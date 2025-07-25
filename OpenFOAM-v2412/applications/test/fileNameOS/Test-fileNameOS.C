/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Test-fileNameOS

Description
    Test fileName behaviour, potential OS capabilities etc.

    In the distant future could possibly replace parts with C++ filesystem

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "OSspecific.H"
#include "Switch.H"

#include <cstdlib>
#include <iostream>


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void testDirname(const std::string& rawInput)
{
    fileName input(fileName::validate(rawInput));

    Info<< nl
        << "input:   " << rawInput << nl
        << "fileName:" << input << nl
        << "   path:" << input.path()
        << "   name:\"" << input.name() << '"'
        << "   ext:\"" << input.ext()  << '"'
        << "   components: " << flatOutput(input.components()) << nl;

    if (rawInput.size() != input.size())
    {
        Info<< "   This would be Fatal with debug > 1" << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption("no-space", "allowSpaceInFileName = false");
    argList::addBoolOption("with-space", "set allowSpaceInFileName = true");

    #include "setRootCase.H"

    if (args.found("with-space"))
    {
        fileName::allowSpaceInFileName = true;
    }

    if (args.found("no-space"))
    {
        fileName::allowSpaceInFileName = false;

    }


    Info<<"fileName with spaces? : "
        << Switch(bool(fileName::allowSpaceInFileName)) << nl << nl;


    {
        testDirname("/abc");
        testDirname("/abc/with space/name");
        testDirname("/abc/with space/more space");
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
