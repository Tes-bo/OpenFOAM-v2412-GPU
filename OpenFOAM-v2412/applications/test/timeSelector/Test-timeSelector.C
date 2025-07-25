/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Test TimePaths and timeSelector

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "Time.H"
#include "timeSelector.H"

using namespace Foam;

bool print(const instantList& instants)
{
    if (instants.empty())
    {
        Info<<"    none" << nl << nl;
        return false;
    }

    Info <<"(" << nl;

    for (const instant& t : instants)
    {
        Info<<"    " << t << nl;
    }

    Info<<")" << nl << nl;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::addNote("Test timeSelector and TimePaths");

    timeSelector::addOptions(true, true);
    argList::noFunctionObjects();

    argList::addOption("relative", "PATH", "Test relativePath");

    #include "setRootCase.H"
    #include "createTime.H"

    Pout<< "Time" << nl
        << "rootPath:   " << runTime.rootPath() << nl
        << "path:       " << runTime.path() << nl
        << "globalCase: " << runTime.globalCaseName() << nl
        << "globalPath: " << runTime.globalPath() << nl
        << nl;

    if (args.found("relative"))
    {
        Pout<< "input path: " << args["relative"] << nl
            << "relative  : " << runTime.relativePath(args["relative"], true)
            << nl
            << nl;
    }

    autoPtr<TimePaths> timePaths;

    instantList times(TimePaths(args).times());

    Info<<"Available times" << nl;

    if (print(times))
    {
        times = timeSelector::select(times, args);

        Info<< "Selected times" << nl;
        print(times);
    }

    return 0;
}

// ************************************************************************* //
