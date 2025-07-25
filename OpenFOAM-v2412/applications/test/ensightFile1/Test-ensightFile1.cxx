/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Test-ensightFile

Description
    check cleanup of ensight file and variable names

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "ensightFile.H"
#include "ensightGeoFile.H"
#include "Switch.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addBoolOption("ascii", "open as ascii instead of binary");
    argList::addBoolOption("binary", "(default)");
    argList::addBoolOption("clear", "force clear of time-steps");
    argList::addBoolOption("no-end", "skip use of endTimeStep");
    argList::addBoolOption("append", "open in append mode");
    argList::addOption("geom", "geometry file");
    argList::addOption("field", "field file");

    #include "setRootCase.H"

    const bool with_ascii = args.found("ascii") && !args.found("binary");
    // const bool with_binary = args.found("binary");
    const bool with_append = args.found("append");
    const bool with_clear = args.found("clear");
    const bool without_end = args.found("no-end");

    const IOstreamOption::streamFormat fmt =
    (
        with_ascii
      ? IOstreamOption::ASCII
      : IOstreamOption::BINARY
    );

    const IOstreamOption::appendType append =
    (
        with_append
      ? IOstreamOption::APPEND_ATE
      : IOstreamOption::NO_APPEND
    );


    fileName file;
    if (args.readIfPresent("geom", file))
    {
        Info<< "Open " << file << " as geometry "
            << " format:" << (with_ascii ? "ASCII" : "BINARY")
            << " append:" << Switch::name(with_append) << nl;

        ensightGeoFile ensFile(append, file, fmt);

        if (append)
        {
            ensFile.beginTimeStep();

            // At the moment need to pair begin/end time-step calls
            if (!without_end)
            {
                ensFile.endTimeStep();
            }
        }

        if (with_clear)
        {
            ensFile.clearTimeSteps();
        }
    }

    if (args.readIfPresent("field", file))
    {
        Info<< "Open " << file << " as field"
            << " format:" << (with_ascii ? "ASCII" : "BINARY")
            << " append:" << Switch::name(with_append) << nl;

        ensightFile ensFile(append, file, fmt);

        if (append)
        {
            ensFile.beginTimeStep();

            // At the moment need to pair begin/end time-step calls
            if (!without_end)
            {
                ensFile.endTimeStep();
            }
        }

        if (with_clear)
        {
            ensFile.clearTimeSteps();
        }
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
