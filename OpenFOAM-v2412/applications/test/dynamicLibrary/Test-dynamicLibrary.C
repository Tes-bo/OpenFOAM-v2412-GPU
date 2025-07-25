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
    Test-dynamicLibrary

Description
    Test loading/unloading of libraries

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "DynamicList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Low-level test of library load/unload");

    profiling::disable(); // No profiling output
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::removeOption("case");
    argList::addBoolOption("no-close", "Skip dlclose");
    argList::addBoolOption("quiet", "Disable verbosity");

    argList::addArgument("lib...");
    argList::noMandatoryArgs();  // Arguments are optional

    argList args(argc, argv, false, true);

    const bool noClose = args.found("no-close");
    const bool verbose = !args.found("quiet");

    //- Pointers to the loaded libraries
    DynamicList<void*> libPtrs_;

    //- Names of loaded libraries, or of libraries to be loaded
    DynamicList<fileName> libNames_;

    label nbad = 0;
    wordHashSet loaded;

    for (int argi = 1; argi < args.size(); ++argi)
    {
        const auto libName = args.get<fileName>(argi);

        if (libName.empty())
        {
            continue;
        }

        void* ptr = Foam::dlOpen(libName, false);

        if (!ptr)
        {
            ++nbad;
        }
        else
        {
            libPtrs_.push_back(ptr);
            libNames_.push_back(libName);

            if (verbose)
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

    if (!noClose)
    {
        forAllReverse(libPtrs_, i)
        {
            void* ptr = libPtrs_[i];

            if (ptr == nullptr)
            {
                libNames_[i].clear();
                continue;
            }

            const bool ok = Foam::dlClose(ptr);

            if (verbose)
            {
                if (ok)
                {
                    InfoErr << "Closed ";
                }
                else
                {
                    InfoErr << "Failed closing ";
                }

                InfoErr
                    << libNames_[i]
                    << " with handle " << Foam::name(ptr) << nl;
            }
        }
    }

    return 0;
}


// ************************************************************************* //
