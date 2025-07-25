/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Test-plotFunction1

Description
    Plot scalar Function1 entries

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fieldTypes.H"
#include "Function1.H"
#include "PtrList.H"
#include "Fstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::setAdvanced("case");  // Hide -case : has no meaning here

    argList::addOption("begin", "scalar", "The start time (default: 0)");
    argList::addOption("end", "scalar", "The end time (default: 1)");
    argList::addOption("incr", "scalar", "The time increment (default: 0.1)");
    argList::addOption("timeBase", "scalar", "The time base (default: 1)");
    argList::addNote
    (
        "Read scalar functions from each file and produces"
        " time/value output for each"
    );

    argList::noMandatoryArgs();
    argList::addArgument("file1");
    argList::addArgument("...");
    argList::addArgument("fileN");

    #include "setRootCase.H"

    const scalar begTime = args.getOrDefault<scalar>("begin", 0);
    const scalar endTime = args.getOrDefault<scalar>("end", 1);
    const scalar incrTime = args.getOrDefault<scalar>("incr", 0.1);
    const scalar timeBase = args.getOrDefault<scalar>("timeBase", 1);

    Info.stream().precision(10);

    for (label argi=1; argi < args.size(); ++argi)
    {
        IFstream is(args.get<fileName>(argi));

        dictionary dict(is);

        for (const entry& dEntry : dict)
        {
            autoPtr<Function1<scalar>> funPtr =
                Function1<scalar>::New(dEntry.keyword(), dict);

            auto& fun = *funPtr;

            InfoErr<< nl;
            fun.writeData(InfoErr);
            InfoErr<< nl;

            Info<< nl << "# " << fun.type() << nl;

            for (scalar t = begTime; t < endTime; t += incrTime)
            {
                Info<< t << tab << fun.value(t*timeBase) << nl;
            }
            Info<< nl;
        }
    }

    return 0;
}


// ************************************************************************* //
