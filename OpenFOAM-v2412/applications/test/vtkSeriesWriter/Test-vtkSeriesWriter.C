/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Test-vtkSeriesWriter

Description
    Basic functionality tests for vtk::seriesWriter

\*---------------------------------------------------------------------------*/

#include "foamVtkSeriesWriter.H"
#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("sort", "Sort value / name");
    argList::addBoolOption("check", "Check for existence of files");
    argList::addOption("time", "value", "Filter based on given time");
    argList::addOption("scan", "series", "Scan directory to create series");

    argList args(argc, argv, false, true);

    const scalar currTime = args.getOrDefault<scalar>("time", GREAT);

    Info<< "Using currTime = " << currTime << nl
        << "when loading " << (args.size()-1) << " files" << nl << nl;

    for (label argi=1; argi < args.size(); ++argi)
    {
        const auto input = args.get<fileName>(argi);

        Info << "load from " << input << nl;

        vtk::seriesWriter writer;
        writer.load(input);

        writer.print(Info);
        Info<< nl << nl;

        if (writer.removeNewer(currTime))
        {
            Info<< "removed entries with time >= " << currTime << nl;
            writer.print(Info);
            Info<< nl << nl;
        }

        if (args.found("sort"))
        {
            writer.sort();

            Info<< "sorted" << nl;
            writer.print(Info);
            Info<< nl << nl;
        }

        if (args.found("check"))
        {
            writer.load(input, true);

            Info<< "reload, checking the existence of files" << nl;
            writer.print(Info);
            Info<< nl << nl;
        }

        if (writer.empty())
        {
            Info<< "No entries" << nl;
        }
        else
        {
            Info<< writer.size() << " entries" << nl;
        }
    }

    if (args.found("scan"))
    {
        vtk::seriesWriter writer;

        writer.scan(args.get<fileName>("scan"));

        Info<< "scanned for files" << nl;
        writer.print(Info);
        Info<< nl << nl;
    }


    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
