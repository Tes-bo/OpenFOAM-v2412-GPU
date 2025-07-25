/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    profilingSummary

Group
    grpMiscUtilities

Description
    Collects information from profiling files in the processor
    sub-directories and summarizes the number of calls and time spent as
    max/avg/min values. If the values are identical for all processes,
    only a single value is written.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "polyMesh.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "OFstream.H"
#include "argList.H"
#include "stringOps.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "functionObject.H"

using namespace Foam;

// The name of the sub-dictionary entry for profiling fileName:
static const word profilingFileName("profiling");

// The name of the sub-dictionary entry for profiling:
static const word blockNameProfiling("profiling");

// The name of the sub-dictionary entry for profiling and tags of entries
// that will be processed to determine (max,avg,min) values
const HashTable<wordList> processing
{
    { "profiling", { "calls", "totalTime", "childTime", "maxMem" } },
    { "memInfo", { "size", "free" } },
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Collect profiling information from processor directories and"
        " summarize time spent and number of calls as (max avg min) values."
    );

    timeSelector::addOptions(true, true);  // constant(true), zero(true)
    argList::noParallel();
    argList::noFunctionObjects();  // Never use function objects

    // Note that this should work without problems when profiling is active,
    // since we don't trigger it anywhere

    #include "setRootCase.H"
    #include "createTime.H"

    // Determine the processor count
    const label nProcs = fileHandler().nProcs(args.path());

    // Create the processor databases
    PtrList<Time> databases(nProcs);

    forAll(databases, proci)
    {
        databases.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/("processor" + Foam::name(proci)),
                args.allowFunctionObjects(),
                args.allowLibs()
            )
        );
    }

    if (!nProcs)
    {
        FatalErrorInFunction
            << "No processor* directories found"
            << exit(FatalError);
    }


    // Use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    if (timeDirs.empty())
    {
        WarningInFunction
            << "No times selected" << nl << endl;
        return 1;
    }

    // ----------------------------------------------------------------------

    // Processor local profiling information
    List<dictionary> profiles(nProcs);

    // Loop over all times
    forAll(timeDirs, timei)
    {
        // Set time for global database
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        // Name/location for the output summary
        const fileName outputName
        {
            functionObject::outputPrefix,
            "profiling",
            runTime.timeName(),
            profilingFileName
        };


        label nDict = 0;

        // Set time for all databases
        forAll(databases, proci)
        {
            profiles[proci].clear();
            databases[proci].setTime(timeDirs[timei], timei);

            // Look for "uniform/profiling" in each processor directory
            IOobjectList objects
            (
                databases[proci].time(),
                databases[proci].timeName(),
                "uniform"
            );

            const IOobject* ioptr = objects.findObject(profilingFileName);
            if (ioptr)
            {
                IOdictionary dict(*ioptr);

                // Full copy
                profiles[proci] = dict;

                // Assumed to be good if it has 'profiling' sub-dict

                const dictionary* ptr = dict.findDict(blockNameProfiling);
                if (ptr)
                {
                    ++nDict;
                }
            }

            if (nDict < proci)
            {
                break;
            }
        }

        if (nDict != nProcs)
        {
            Info<< "found " << nDict << "/" << nProcs
                << " profiling files" << nl << endl;
            continue;
        }


        // Information seems to be there for all processors
        // can do a summary

        IOdictionary summary
        (
            IOobject
            (
                runTime.path()/outputName,
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER,
                true   // global-like
            )
        );

        summary.note() =
        (
            "summarized (max avg min) values from "
          + Foam::name(nProcs) + " processors"
        );


        // Accumulator for each tag
        HashTable<DynamicList<scalar>> stats;

        // Use first as 'master' to decide what others have
        forAllConstIters(profiles.first(), mainIter)
        {
            const entry& mainEntry = mainIter();

            // level1: eg, profiling {} or memInfo {}
            const word& level1Name = mainEntry.keyword();

            if
            (
                !processing.found(level1Name)
             || !mainEntry.isDict()
             || mainEntry.dict().empty()
            )
            {
                continue;  // Only process known types
            }

            const wordList& tags = processing[level1Name];

            const dictionary& level1Dict = mainEntry.dict();

            // We need to handle sub-dicts with other dicts
            //     Eg, trigger0 { .. } trigger1 { .. }
            //
            // and ones with primitives
            //     Eg, size xx; free yy;

            // Decide based on the first entry:

            // level2: eg, profiling { trigger0 { } }
            // or simply itself it contains primitives only

            wordList level2Names;

            const bool hasDictEntries
                = mainEntry.dict().first()->isDict();

            if (hasDictEntries)
            {
                level2Names =
                    mainEntry.dict().sortedToc(stringOps::natural_sort());
            }
            else
            {
                level2Names = {level1Name};
            }

            summary.set(level1Name, dictionary());

            dictionary& outputDict = summary.subDict(level1Name);

            for (const word& level2Name : level2Names)
            {
                // Presize everything
                stats.clear();
                for (const word& tag : tags)
                {
                    stats(tag).reserve(nProcs);
                }

                label nEntry = 0;

                for (const dictionary& procDict : profiles)
                {
                    const dictionary* inDictPtr = procDict.findDict(level1Name);

                    if (inDictPtr && hasDictEntries)
                    {
                        // Descend to the next level as required
                        inDictPtr = inDictPtr->findDict(level2Name);
                    }

                    if (!inDictPtr)
                    {
                        break;
                    }

                    ++nEntry;

                    for (const word& tag : tags)
                    {
                        scalar val;

                        if
                        (
                            inDictPtr->readIfPresent(tag, val, keyType::LITERAL)
                        )
                        {
                            stats(tag).append(val);
                        }
                    }
                }

                if (nEntry != nProcs)
                {
                    continue;
                }

                dictionary* outDictPtr = nullptr;

                // Make a full copy of this entry prior to editing it
                if (hasDictEntries)
                {
                    outputDict.add(level2Name, level1Dict.subDict(level2Name));
                    outDictPtr = outputDict.findDict(level2Name);
                }
                else
                {
                    // merge into existing (empty) dictionary
                    summary.add(level1Name, level1Dict, true);
                    outDictPtr = &outputDict;
                }

                dictionary& outSubDict = *outDictPtr;

                // Remove trailing 'processor0' from any descriptions
                // (looks nicer)
                {
                    const word key("description");
                    string val;

                    if (outSubDict.readIfPresent(key, val))
                    {
                        if (val.removeEnd("processor0"))
                        {
                            outSubDict.set(key, val);
                        }
                    }
                }

                // Process each tag (calls, time etc)
                for (const word& tag : tags)
                {
                    DynamicList<scalar>& lst = stats(tag);

                    if (lst.size() == nProcs)
                    {
                        sort(lst);
                        const scalar avg = sum(lst) / nProcs;

                        if (lst.first() != lst.last())
                        {
                            outSubDict.set
                            (
                                tag,
                                FixedList<scalar, 3>
                                {
                                    lst.last(), avg, lst.first()
                                }
                            );
                        }
                    }
                }
            }
        }


        // Now write the summary
        {
            mkDir(summary.path());

            OFstream os(summary.objectPath());

            summary.writeHeader(os);
            summary.writeData(os);
            IOobject::writeEndDivider(os);

            Info<< "Wrote to " << outputName << nl << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
