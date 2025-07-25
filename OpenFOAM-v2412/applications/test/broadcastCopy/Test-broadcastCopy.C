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

Description
    Test file/directory broadcasting

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"
#include "fileOperation.H"
#include "Pstream.H"
#include "Switch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noFunctionObjects();
    argList::noCheckProcessorDirectories();

    argList::addNote("Test broadcast file via MPI");

    argList::addArgument("srcFile");
    argList::addBoolOption("even", "Broadcast to even directories only");
    argList::addBoolOption("relative", "Copy relative to output dir");

    #include "setRootCase.H"

    const auto srcFile = args.get<fileName>(1);

    // const auto dstFile = args.get<fileName>(2);
    fileName dstFile("proc" + Foam::name(UPstream::myProcNo()));


    const bool writeOnProc =
    (
        !args.found("even") || 0 == (UPstream::myProcNo() % 2)
    );

    if (args.found("relative"))
    {
        // if (writeOnProc)
        // {
        //     Foam::mkDir(dstFile);
        // }
    }
    else
    {
        dstFile /= srcFile + ".copy";
    }

    Pout<< "writing: " << writeOnProc << " : " << dstFile << endl;

    const auto& fp = fileHandler();

    fp.broadcastCopy
    (
        UPstream::worldComm,
        writeOnProc,
        srcFile,
        dstFile
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
