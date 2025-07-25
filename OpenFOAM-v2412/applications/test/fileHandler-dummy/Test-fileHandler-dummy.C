/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-fileHandler-dummy

Description
    Simple test of dummy fileOperation

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "fileOperation.H"
#include "Switch.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();

    argList::addBoolOption
    (
        "force",
        "Force use of dummy handler (and provoke NotImplemented)"
    );

    #include "setRootCase.H"

    const bool optForce = args.found("force");

    const auto& dummy = fileOperation::null();

    Info<< "default handler: " << fileHandler() << endl;
    Info<< "dummy handler: " << dummy() << endl;

    Switch hasFile(Switch::INVALID);

    if (optForce || (dummy && dummy().good()))
    {
        hasFile = dummy().isFile("foo");
    }

    Info<< "check file: " << hasFile << endl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
