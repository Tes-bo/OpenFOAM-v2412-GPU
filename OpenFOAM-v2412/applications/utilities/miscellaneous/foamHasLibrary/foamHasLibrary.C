/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    foamHasLibrary

Group
    grpMiscUtilities

Description
    Test if given libraries can be loaded.

Usage
    \b foamHasLibrary [OPTION] lib...

    Options:
      - \par -or
        Success if any of the libraries can be loaded.
        Does not short-circuit.

      - \par -detail
        Additional detail (meaning may change).

      - \par -verbose
        Additional verbosity

Note
    No normal output.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Test if given libraries can be loaded");

    profiling::disable(); // No profiling output
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::removeOption("case");
    argList::addBoolOption
    (
        "or",
        "Success if any of the libraries can be loaded\n"
        "(does not short-circuit)"
    );
    argList::addBoolOption
    (
        "detail",
        "Additional detail"
    );
    argList::addVerboseOption();

    argList::addArgument("lib...");
    argList::noMandatoryArgs();  // Arguments are optional

    argList args(argc, argv, false, true);

    // Force dlOpen of FOAM_DLOPEN_LIBS (principally for Windows applications)
    #include "foamDlOpenLibs.H"

    const bool testOr = args.found("or");
    const bool detail = args.found("detail");

    label ngood = 0;
    label nbad = 0;

    dlLibraryTable& libs = args.libs();

    wordHashSet loaded;

    for (int argi = 1; argi < args.size(); ++argi)
    {
        const auto libName = args.get<fileName>(argi);  // with validate

        if (libName.empty())
        {
            continue;
        }

        // InfoErr << "Check " << libName << nl;

        // Could have libs.findLibrary(...)
        // if we really expect many duplicates

        const void* ptr = libs.open(libName, false);

        if (!ptr)
        {
            ++nbad;
        }
        else
        {
            ++ngood;

            if (args.verbose())
            {
                const word addr(Foam::name(ptr));

                if (loaded.insert(addr))
                {
                    InfoErr << "Can load " << libName << nl;
                }
                else
                {
                    InfoErr << "Already loaded " << libName << nl;
                }
            }
        }
    }

    if (detail)
    {
        InfoErr << libs.info();
    }

    return (nbad == 0 || (testOr && ngood > 0)) ? 0 : 1;
}


// ************************************************************************* //
