/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    foamListTimes

Group
    grpPostProcessingUtilities

Description
    List times using the timeSelector, or use to remove selected time
    directories.

Usage
    \b foamListTimes [OPTION]

    Options:
      - \par -processor
        Times from processor0/ directory

      - \par -rm
        Remove selected time directories

      - \par -verbose
        Report progress during removal

Note
    The OpenFOAM banner information is suppressed so that the output can be
    piped into another command.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "autoPtr.H"
#include "profiling.H"
#include "timeSelector.H"
#include "TimePaths.H"
#include "ListOps.H"
#include "stringOps.H"

using namespace Foam;

// Many ways to name processor directories
//
// Uncollated       | "processor0", "processor1" ...
// Collated         | "processors<N>"
// Host collated    | "processors<N>_<low>-<high>"

const regExp matcher("processors?[0-9]+(_[0-9]+-[0-9]+)?");

bool isProcessorDir(const string& dir)
{
    return (dir.starts_with("processor") && matcher.match(dir));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "List times using the timeSelector,"
        " or use to remove selected time directories"
    );
    timeSelector::addOptions(true, true);  // constant(true), zero(true)
    argList::noBanner();
    argList::noParallel();
    argList::noJobInfo();
    argList::noFunctionObjects();  // Never use function objects
    argList::addBoolOption
    (
        "processor",
        "List times from processor0/ directory"
    );
    argList::addBoolOption
    (
        "rm",
        "Remove selected time directories"
    );
    argList::addVerboseOption
    (
        "Report progress of -rm option"
    );
    profiling::disable(); // Disable profiling (and its output)

    #include "setRootCase.H"

    const bool removeFiles(args.found("rm"));
    bool verbose(args.verbose());


    // Get times list from the master processor and subset based on
    // command-line options

    label nProcs = 0;
    autoPtr<TimePaths> timePaths;

    if (args.found("processor"))
    {
        // Determine the processor count
        nProcs = fileHandler().nProcs(args.path());

        if (!nProcs)
        {
            FatalErrorInFunction
                << "No processor* directories found"
                << exit(FatalError);
        }

        // Obtain time directory names from "processor0/" only
        timePaths = autoPtr<TimePaths>::New
        (
            args.rootPath(),
            args.caseName()/"processor0"
        );
    }
    else
    {
        timePaths = autoPtr<TimePaths>::New
        (
            args.rootPath(),
            args.caseName()
        );
    }


    const instantList timeDirs(timeSelector::select(timePaths->times(), args));

    const label nTimes = timeDirs.size();

    if (removeFiles)
    {
        if (nProcs)
        {
            fileNameList procDirs
            (
                Foam::readDir
                (
                    args.path(),
                    fileName::DIRECTORY,
                    false,  // No gzip anyhow
                    false   // Do not follow linkts
                )
            );

            inplaceSubsetList(procDirs, isProcessorDir);

            // Perhaps not needed
            /// Foam::sort(procDirs, stringOps::natural_sort());

            if (verbose)
            {
                InfoErr
                    << "Removing " << nTimes
                    << " times in " << procDirs.size()
                    << " processor directories" << endl;
            }

            // No processor directories? - silence verbosity
            if (procDirs.empty())
            {
                verbose = false;
            }

            forAllReverse(timeDirs, timei)
            {
                const word& timeName = timeDirs[timei].name();

                if (verbose)
                {
                    InfoErr
                        << "    rm " << timeName
                        << " [" << (nTimes - timei) << '/' << nTimes << ']'
                        << endl;
                }

                for (const fileName& procDir : procDirs)
                {
                    rmDir(args.path()/procDir/timeName, true);
                }
            }
        }
        else
        {
            if (verbose)
            {
                InfoErr
                    << "Removing " << nTimes
                    << " time directories" << endl;
            }

            forAllReverse(timeDirs, timei)
            {
                const word& timeName = timeDirs[timei].name();

                if (verbose)
                {
                    InfoErr
                        << "    rm " << timeName
                        << " [" << (nTimes - timei) << '/' << nTimes << ']'
                        << endl;
                }

                rmDir(args.path()/timeName, true);
            }
        }
    }
    else
    {
        // List times: one per line
        for (const instant& t : timeDirs)
        {
            Info<< t.name() << nl;
        }
        Info<< flush;
    }

    return 0;
}


// ************************************************************************* //
