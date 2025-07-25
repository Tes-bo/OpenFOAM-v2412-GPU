/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    foamRestoreFields

Group
    grpMiscUtilities

Description
    Adjust (restore) field names by removing the ending.
    The fields are selected automatically or can be specified as optional
    command arguments.

    The operation 'mean' renames files ending with 'Mean' and makes
    a backup of existing names, using the '.orig' ending.

    The operation 'orig' renames files ending with '.orig'.

Usage
    \b foamRestoreFields [OPTION]

    Options:
      - \par -method mean | orig
        The renaming method.

      - \par -processor
        Use processor directories, taking information from processor0/

      - \par -dry-run
        Test without actually moving/renaming files.

      - \par -verbose
        Additional verbosity.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "autoPtr.H"
#include "profiling.H"
#include "timeSelector.H"
#include "Enum.H"
#include "TimePaths.H"
#include "ListOps.H"
#include "stringOps.H"
#include "regionProperties.H"
#include "polyMesh.H"
#include "Time.H"

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


//- The known and support types of operations
enum restoreMethod
{
    MEAN,
    ORIG
};


static const Enum<restoreMethod> methodNames
{
    { restoreMethod::MEAN, "mean" },
    { restoreMethod::ORIG, "orig" },
};


static const Enum<restoreMethod> methodEndings
{
    { restoreMethod::MEAN, "Mean" },
    { restoreMethod::ORIG, ".orig" },
};


// Files in given directory at time instant
inline wordList getFiles(const fileName& dir, const word& instance)
{
    return ListOps::create<word>
    (
        Foam::readDir(dir/instance, fileName::FILE),
        nameOp<fileName>()
    );
}


// Command-line options: -dry-run, -verbose
bool dryrun = false, verbose = false;


// Use predefined method to walk the directory and rename the files.
//
// If no target names are specified, the existing files are scanned for
// candidates.
label restoreFields
(
    const restoreMethod method,
    const fileName& dirName,
    const wordHashSet& existingFiles,
    const wordList& targetNames
)
{
    // The file ending to search for.
    const word ending(methodEndings[method]);

    // The backup ending for existing (if any)
    word bak;

    switch (method)
    {
        case restoreMethod::MEAN:
            bak = methodEndings[restoreMethod::ORIG];
            break;

        default:
            break;
    }

    wordHashSet targets(targetNames);

    if (targets.empty())
    {
        // No target names specified - scan existing files for candidates.

        for (word f : existingFiles)    // Operate on a copy
        {
            // Eg, check for "UMean" and save as "U"
            if (f.removeEnd(ending) && f.size())
            {
                targets.insert(f);
            }
        }
    }

    if (verbose)
    {
        Info<< "directory " << dirName.name() << nl;
    }

    // Count of files moved, including backups
    label count = 0;

    for (const word& dst : targets)
    {
        const word src(dst + ending);

        if (!existingFiles.found(src))
        {
            continue;
        }

        if (bak.size() && existingFiles.found(dst))
        {
            if (dryrun || Foam::mv(dirName/dst, dirName/dst + bak))
            {
                Info<< "    mv  " << dst << "  " << word(dst + bak) << nl;
                ++count;
            }
        }

        if (dryrun || Foam::mv(dirName/src, dirName/dst))
        {
            Info<< "    mv  " << src << "  " << dst << nl;
            ++count;
        }
    }

    return count;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Restore field names by removing the ending. Fields are selected"
        " automatically or can be specified as optional command arguments"
    );

    profiling::disable(); // Disable profiling (and its output)
    argList::noJobInfo();
    argList::noFunctionObjects();  // Never use function objects
    argList::addOption
    (
        "method",
        "name",
        "The restore method (mean|orig) [MANDATORY]. "
        "With <mean> renames files ending with 'Mean' "
        "(with backup of existing as '.orig'). "
        "With <orig> renames files ending with '.orig'"
    );
    argList::addBoolOption
    (
        "processor",
        "In serial mode use times from processor0/ directory, but operate on "
        "processor\\d+ directories"
    );
    argList::addDryRunOption
    (
        "Report action without moving/renaming"
    );
    argList::addVerboseOption();

    // Arguments are optional (non-mandatory)
    argList::noMandatoryArgs();
    argList::addArgument("fieldName ... fieldName");

    timeSelector::addOptions(true, true);  // constant(true), zero(true)

    #include "addAllRegionOptions.H"
    #include "setRootCase.H"

    wordList regionNames0;
    {
        // Dummy time just for the database to read regionProperties

        autoPtr<Time> dummyTimePtr(Time::New(args));

        const auto& runTime = *dummyTimePtr;

        // Handle -allRegions, -regions, -region
        #include "getAllRegionOptions.H"

        regionNames0 = std::move(regionNames);
    }

    dryrun = args.dryRun();
    verbose = args.verbose();


    // Construct time
    // ~~~~~~~~~~~~~~

    restoreMethod method = restoreMethod::ORIG;
    {
        word methodName;

        if
        (
            args.readIfPresent("method", methodName)
         && methodNames.found(methodName)
        )
        {
            method = methodNames[methodName];
        }
        else
        {
            Info<< "Unspecified or unknown method name" << nl
                << "Valid methods: "
                << flatOutput(methodNames.sortedToc()) << nl
                << "... stopping" << nl << nl;
                return 1;
        }
    }

    // Optional base or target field names (eg, 'U', 'T' etc)
    wordList targetNames;
    if (args.size() > 1)
    {
        targetNames.resize(args.size()-1);
        wordHashSet uniq;

        for (label argi=1; argi < args.size(); ++argi)
        {
            if (uniq.insert(args[argi]))
            {
                targetNames[uniq.size()-1] = args[argi];
            }
        }

        targetNames.resize(uniq.size());

        if (verbose)
        {
            Info<< nl
                << "using method=" << methodNames[method] << nl
                << "with fields " << flatOutput(targetNames) << nl;
        }
    }
    else if (verbose)
    {
        Info<< nl
            << "using method=" << methodNames[method] << nl
            << "autodetect fields" << nl;
    }


    // Get times list from the master processor and subset based on
    // command-line options

    label nProcs = 0;
    autoPtr<TimePaths> timePaths;

    if (args.found("processor") && !Pstream::parRun())
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

    fileNameList procDirs;
    label leadProcIdx = -1;

    if (timeDirs.empty())
    {
        Info<< "No times selected" << nl;
    }
    else if (nProcs)
    {
        procDirs =
            Foam::readDir
            (
                args.path(),
                fileName::DIRECTORY,
                false,  // No gzip anyhow
                false   // Do not follow linkts
            );

        inplaceSubsetList(procDirs, isProcessorDir);

        // Perhaps not needed
        Foam::sort(procDirs, stringOps::natural_sort());

        // Decide who will be the "leading" processor for obtaining names
        // - processor0
        // - processors<N>
        // - processors<N>_0-<high>

        // Uncollated
        leadProcIdx = procDirs.find("processor0");

        if (!procDirs.empty())
        {
            if (leadProcIdx < 0)
            {
                // Collated
                leadProcIdx = procDirs.find("processors" + Foam::name(nProcs));
            }

            if (leadProcIdx < 0)
            {
                // Host-collated
                const std::string prefix
                (
                    "processors" + Foam::name(nProcs) + "_0-"
                );

                forAll(procDirs, idx)
                {
                    if (procDirs[idx].starts_with(prefix))
                    {
                        leadProcIdx = idx;
                        break;
                    }
                }
            }

            // Just default to anything (safety)
            if (leadProcIdx < 0)
            {
                leadProcIdx = 0;
            }
        }
    }


    for (const instant& t : timeDirs)
    {
        const word& timeName = t.name();

        Info<< "\nTime = " << timeName << nl;

        for (const word& regionName : regionNames0)
        {
            const word& regionDir = polyMesh::regionName(regionName);

            if (regionNames0.size() > 1)
            {
                Info<< "region = " << regionName << nl;
            }

            label count = 0;
            wordList files;

            if (nProcs)
            {
                if (leadProcIdx >= 0)
                {
                    files =
                        getFiles
                        (
                            args.path()/procDirs[leadProcIdx],
                            timeName/regionDir
                        );
                }

                for (const fileName& procDir : procDirs)
                {
                    count += restoreFields
                    (
                        method,
                        args.path()/procDir/timeName/regionDir,
                        wordHashSet(files),
                        targetNames
                    );
                }
            }
            else
            {
                if (Pstream::master())
                {
                    files = getFiles(args.path(), timeName/regionDir);
                }
                Pstream::broadcast(files);

                count += restoreFields
                (
                    method,
                    args.path()/timeName/regionDir,
                    wordHashSet(files),
                    targetNames
                );
            }

            if (dryrun)
            {
                Info<< "dry-run: ";
            }
            Info<< "moved " << count << " files" << nl;
        }
    }
    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
