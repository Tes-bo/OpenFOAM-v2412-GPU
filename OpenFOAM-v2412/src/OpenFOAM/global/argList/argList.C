/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"
#include "Switch.H"
#include "clock.H"
#include "dictionary.H"
#include "IOobject.H"
#include "JobInfo.H"
#include "labelList.H"
#include "IOobject.H"
#include "dynamicCode.H"
#include "simpleObjectRegistry.H"
#include "sigFpe.H"
#include "sigInt.H"
#include "sigQuit.H"
#include "sigSegv.H"
#include "foamVersion.H"
#include "stringOps.H"
#include "CStringList.H"
#include "stringListOps.H"
#include "fileOperation.H"
#include "fileOperationInitialise.H"

#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::argList::argsMandatory_ = true;
bool Foam::argList::checkProcessorDirectories_ = true;
bool Foam::argList::parallelThreads_ = false;

Foam::SLList<Foam::string>    Foam::argList::validArgs;
Foam::HashSet<Foam::string>   Foam::argList::advancedOptions;
Foam::HashTable<Foam::string> Foam::argList::validOptions;
Foam::HashTable<Foam::string> Foam::argList::validParOptions;
Foam::HashTable<Foam::string, Foam::label, Foam::Hash<Foam::label>>
    Foam::argList::argUsage;
Foam::HashTable<Foam::string> Foam::argList::optionUsage;
Foam::HashTable<std::pair<Foam::word,int>> Foam::argList::validOptionsCompat;
Foam::HashTable<std::pair<bool,int>> Foam::argList::ignoreOptionsCompat;
Foam::SLList<Foam::string>    Foam::argList::notes;
Foam::word Foam::argList::postProcessOptionName("postProcess");

std::string::size_type Foam::argList::usageMin = 20;
std::string::size_type Foam::argList::usageMax = 80;

Foam::argList::initValidTables::initValidTables()
{
    argList::addOption
    (
        "case",
        "dir",
        "Case directory (instead of current directory)"
    );
    argList::addOption
    (
        "lib",
        "name",
        "Additional library or library list to load.\n"
        "[Can be used multiple times]",
        true  // advanced option
    );
    argList::addBoolOption
    (
        "no-libs",
        "Disable use of the controlDict 'libs' entry",
        true  // advanced option
    );

    argList::addOption
    (
        "debug-switch",
        "name=val",
        "Set named DebugSwitch (default value: 1).\n"
        "[Can be used multiple times]",
        true  // advanced option
    );

    argList::addOption
    (
        "info-switch",
        "name=val",
        "Set named InfoSwitch (default value: 1).\n"
        "[Can be used multiple times]",
        true  // advanced option
    );

    argList::addOption
    (
        "opt-switch",
        "name=val",
        "Set named OptimisationSwitch (default value: 1).\n"
        "[Can be used multiple times]",
        true  // advanced option
    );

    argList::addBoolOption("parallel", "Run in parallel");
    validParOptions.set("parallel", "");

    argList::addBoolOption
    (
        "mpi-threads",
        "Request use of MPI threads",
        true  //  advanced option
    );

    argList::addOption
    (
        "roots",
        "(dir1 .. dirN)",
        "Subprocess root directories for distributed running",
        true  // advanced option
    );
    validParOptions.set
    (
        "roots",
        "(dir1 .. dirN)"
    );

    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Alternative decomposePar dictionary file"
    );
    argList::addOption
    (
        "hostRoots",
        "((host1 dir1) .. (hostN dirN))",
        "Per-subprocess root directories for distributed running."
        " The host specification can be a regex.",
        true  // advanced option
    );
    validParOptions.set
    (
        "hostRoots",
        "((host1 dir1) .. (hostN dirN))"
    );

    argList::addBoolOption
    (
        "noFunctionObjects",
        "Do not execute function objects",
        true  // advanced option
    );

    argList::addOption
    (
        "fileHandler",
        "handler",
        "Override the file handler type",
        true  // advanced option
    );

    argList::addOption
    (
        "world",
        "name",
        "Name of the local world for parallel communication",
        true  // advanced option
    );
    validParOptions.set
    (
        "world",
        "name"
    );


    // Some standard option aliases (with or without version warnings)
//     argList::addOptionCompat
//     (
//         "noFunctionObjects", {"no-function-objects", 0}
//     );

    Pstream::addValidParOptions(validParOptions);
}

Foam::argList::initValidTables dummyInitValidTables;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Counted per machine name
// Does not include any sorting since we wish to know the ordering according to
// mpi rank.
//
// Always include the master too.
// This provides a better overview of the subscription
static void printHostsSubscription(const UList<string>& hostProcs)
{
    Info<< "Hosts  :\n(" << nl;

    std::string prev = Foam::hostName();
    int count = 1;

    for (const auto& str : hostProcs)
    {
        std::string curr(str.substr(0, str.rfind('.')));

        if (prev != curr)
        {
            if (count)
            {
                // Finish previous
                Info<< "    (" << prev.c_str() << ' ' << count << ')' << nl;
                count = 0;
            }

            prev = std::move(curr);
        }
        ++count;
    }

    if (count)
    {
        // Finished last one
        Info<< "    (" << prev.c_str() << ' ' << count << ')' << nl;
    }

    Info<< ')' << nl;
}


static bool printRootsSubscription
(
    const UList<string>& hostProcs,
    const UList<fileName>& roots
)
{
    if (hostProcs.size() == roots.size()+1)
    {
        // Sort roots according to hostProc
        DynamicList<string> sortedProcs;
        DynamicList<label> sortedRoots;

        forAll(roots, i)
        {
            const fileName& root = roots[i];
            const string& host = hostProcs[i+1];
            const label index = sortedProcs.find(host);

            if (index == -1)
            {
                sortedProcs.push_back(host);
                sortedRoots.push_back(i);
            }
            else if (roots[sortedRoots[index]] != root)
            {
                // Not properly sorted...
                return false;
            }
        }

        Info<< "Roots  :\n(" << nl;
        forAll(sortedProcs, i)
        {
            Info<< "    (" << sortedProcs[i].c_str() << ' '
                << roots[sortedRoots[i]] << ')' << nl;
        }
        Info<< ')' << nl;

        return true;
    }

    return false;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::argList::checkITstream(const ITstream& is, const label index)
{
    const label remaining = is.nRemainingTokens();

    if (remaining)
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Argument " << index << " has "
            << remaining << " excess tokens" << nl << nl;
    }
    else if (!is.size())
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Argument " << index << " had no tokens" << nl << nl;
    }
}


void Foam::argList::checkITstream(const ITstream& is, const word& optName)
{
    const label remaining = is.nRemainingTokens();

    if (remaining)
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Option -" << optName << " has "
            << remaining << " excess tokens" << nl << nl;
    }
    else if (!is.size())
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Option -" << optName << " had no tokens" << nl << nl;
    }
}


void Foam::argList::raiseBadInput(const word& optName) const
{
    // Can use FatalError
    // predicate checks are not used at the earliest stages
    FatalErrorIn(executable())
        << "Option -" << optName << " with invalid input" << nl
        << exit(FatalError);
}


void Foam::argList::addArgument
(
    const string& argName,
    const string& usage
)
{
    validArgs.push_back(argName);

    // The first program argument starts at 1 - obtain index after push_back()

    const label index = validArgs.size();

    if (usage.empty())
    {
        argUsage.erase(index);
    }
    else
    {
        argUsage.set(index, usage);
    }
}


void Foam::argList::addBoolOption
(
    const word& optName,
    const string& usage,
    bool advanced
)
{
    argList::addOption(optName, "", usage, advanced);
}


void Foam::argList::addOption
(
    const word& optName,
    const string& param,
    const string& usage,
    bool advanced
)
{
    validOptions.set(optName, param);
    if (!usage.empty())
    {
        optionUsage.set(optName, usage);
    }
    if (advanced)
    {
        advancedOptions.set(optName);
    }
}


void Foam::argList::setAdvanced(const word& optName, bool advanced)
{
    if (advanced && validOptions.contains(optName))
    {
        advancedOptions.set(optName);
    }
    else
    {
        advancedOptions.erase(optName);
    }
}


void Foam::argList::addOptionCompat
(
    const word& optName,
    std::pair<const char*,int> compat
)
{
    validOptionsCompat.insert
    (
        compat.first,
        std::pair<word,int>(optName, compat.second)
    );
}


void Foam::argList::ignoreOptionCompat
(
    std::pair<const char*,int> compat,
    bool expectArg
)
{
    ignoreOptionsCompat.insert
    (
        compat.first,
        std::pair<bool,int>(expectArg, compat.second)
    );
}


void Foam::argList::addUsage
(
    const word& optName,
    const string& usage
)
{
    if (usage.empty())
    {
        optionUsage.erase(optName);
    }
    else
    {
        optionUsage.set(optName, usage);
    }
}


void Foam::argList::addNote(const string& note)
{
    if (!note.empty())
    {
        notes.push_back(note);
    }
}


void Foam::argList::removeOption(const word& optName)
{
    validOptions.erase(optName);
    optionUsage.erase(optName);
    advancedOptions.erase(optName);
}


void Foam::argList::noMandatoryArgs()
{
    argsMandatory_ = false;
}


bool Foam::argList::argsMandatory()
{
    return argsMandatory_;
}


void Foam::argList::noBanner()
{
    ::Foam::infoDetailLevel = 0;
}


bool Foam::argList::bannerEnabled()
{
    return (::Foam::infoDetailLevel > 0);
}


void Foam::argList::addDryRunOption
(
    const string& usage,
    bool advanced
)
{
    const word optName("dry-run", false);

    argList::addBoolOption(optName, usage, advanced);
    if (!advanced)
    {
        advancedOptions.erase(optName);  // Avoid 'stickiness'
    }
}


void Foam::argList::addVerboseOption
(
    const string& usage,
    bool advanced
)
{
    const word optName("verbose", false);

    if (usage.empty())
    {
        argList::addBoolOption
        (
            optName,
            "Additional verbosity (can be used multiple times)",
            advanced
        );
    }
    else
    {
        argList::addBoolOption(optName, usage, advanced);
    }
    if (!advanced)
    {
        advancedOptions.erase(optName);  // Avoid 'stickiness'
    }
}


void Foam::argList::noFunctionObjects(bool addWithOption)
{
    removeOption("noFunctionObjects");

    // Ignore this bool option without warning
    // - cannot tie to any particular version anyhow
    ignoreOptionCompat({"noFunctionObjects", 0}, false);

    if (addWithOption)
    {
        argList::addBoolOption
        (
            "withFunctionObjects",
            "Execute functionObjects"
            // An advanced option, but seldom used so expose it more
        );
    }
}


void Foam::argList::noJobInfo()
{
    JobInfo::disable();
}


void Foam::argList::noLibs()
{
    argList::addBoolOption
    (
        "no-libs",
        "Disable use of the controlDict libs entry",
        true  // advanced option
    );
}


void Foam::argList::noParallel()
{
    removeOption("parallel");
    removeOption("roots");
    removeOption("decomposeParDict");
    removeOption("hostRoots");
    removeOption("world");
    removeOption("mpi-threads");
    validParOptions.clear();
}


void Foam::argList::parallelThreads_on()
{
    parallelThreads_ = true;
}


void Foam::argList::noCheckProcessorDirectories()
{
    checkProcessorDirectories_ = false;
}


bool Foam::argList::postProcess(int argc, char *argv[])
{
    for (int argi = 1; argi < argc; ++argi)
    {
        const char *optName = argv[argi];

        if (optName[0] == '-')
        {
            ++optName;  // Looks like an option, skip leading '-'

            if (optName == postProcessOptionName)
            {
                return true;
            }
        }
    }

    return false;
}


int Foam::argList::verbose(int argc, char *argv[])
{
    int num = 0;
    for (int argi = 1; argi < argc; ++argi)
    {
        if (strcmp(argv[argi], "-verbose") == 0)
        {
            ++num;
        }
    }
    return num;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::argList::envExecutable()
{
    return Foam::getEnv("FOAM_EXECUTABLE");
}


Foam::fileName Foam::argList::envGlobalPath()
{
    return Foam::getEnv("FOAM_CASE");
}


Foam::fileName Foam::argList::envRelativePath
(
    const fileName& input,
    const bool caseTag
)
{
    if (input.isAbsolute())
    {
        return input.relative(envGlobalPath(), caseTag);
    }

    return input;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::argList::optionCompat(const word& optName)
{
    // NB: optName includes the leading '-' so that the return value
    // can be used directly

    if (!validOptionsCompat.empty())
    {
        const auto fnd = validOptionsCompat.cfind(optName.substr(1));

        if (fnd.good())
        {
            const auto& alt = fnd.val();

            // No error::master() guard - only called on master anyhow
            if (error::warnAboutAge(alt.second))
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Found [v" << alt.second << "] '"
                    << optName << "' instead of '-"
                    << alt.first << "' option"
                    << nl
                    << std::endl;

                error::warnAboutAge("option", alt.second);
            }

            return "-" + alt.first;
        }
    }

    // Nothing found - pass through the original input
    return optName;
}


int Foam::argList::optionIgnore(const word& optName)
{
    // NB: optName is without the leading '-'

    if (!ignoreOptionsCompat.empty())
    {
        const auto fnd = ignoreOptionsCompat.cfind(optName);

        if (fnd.good())
        {
            const auto& alt = fnd.val();

            // Number to skip (including the option itself)
            // '-option ARG' or '-option'
            const int nskip = (alt.first ? 2 : 1);

            // No error::master() guard - only called on master anyhow
            if (error::warnAboutAge(alt.second))
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Ignoring [v" << alt.second << "] '-"
                    << optName << (nskip > 1 ? " ARG" : "")
                    << "' option"
                    << nl
                    << std::endl;

                error::warnAboutAge("option", alt.second);
            }

            return nskip;
        }
    }

    // TBD: could ignore -verbose, -dry-run etc if they are not active...

    return 0; // Do not skip
}


bool Foam::argList::regroupArgv(int& argc, char**& argv)
{
    int nArgs = 1;
    int ignore = 0;
    unsigned depth = 0;
    string group;  // For grouping ( ... ) arguments

    // Note: we rewrite directly into args_
    // and use a second pass to sort out args/options

    args_[0] = fileName(argv[0]);
    for (int argi = 1; argi < argc; ++argi)
    {
        const char *optName = argv[argi];

        if (optName[0] == '(' && optName[1] == '\0')
        {
            // Begin list
            ++depth;
            group += '(';
        }
        else if (optName[0] == ')' && optName[1] == '\0')
        {
            if (depth)
            {
                // End list
                --depth;
                group += ')';
                if (!depth)
                {
                    args_[nArgs++] = group;
                    group.clear();
                }
            }
            else
            {
                // A stray ')' - likely never happens
                args_[nArgs++] = argv[argi];
            }
        }
        else if (depth)
        {
            // Quote each string element
            group += '"';
            group += argv[argi];
            group += '"';
        }
        else if (optName[0] == '-')
        {
            ++optName;  // Looks like an option, skip leading '-'

            if (validOptions.contains(optName))
            {
                // Known option name
                args_[nArgs++] = argv[argi];
            }
            else if ((ignore = optionIgnore(optName)) > 0)
            {
                // Option to be ignored (with/without an argument)
                if (ignore > 1)
                {
                    ++argi;
                }
            }
            else
            {
                // Try alias for the option name
                args_[nArgs++] = optionCompat(argv[argi]);
            }
        }
        else
        {
            args_[nArgs++] = argv[argi];
        }
    }

    if (group.size())
    {
        // Group(s) not closed, but flush anything still pending
        args_[nArgs++] = group;
    }

    args_.resize(nArgs);

    std::string::size_type len = (nArgs-1); // Spaces between args
    for (const auto& s : args_)
    {
        len += s.length();
    }

    // Length needed for regrouped command-line
    commandLine_.reserve(len);

    return nArgs < argc;
}


void Foam::argList::setCasePaths()
{
    fileName caseDir;

    const auto optIter = options_.cfind("case");  // [-case dir] specified?

    if (optIter.good())
    {
        caseDir = fileName::validate(optIter.val());  // includes 'clean'

        if (caseDir.empty() || caseDir == ".")
        {
            // Treat "", "." and "./" as if -case was not specified
            caseDir = cwd();
            options_.erase("case");
        }
        else
        {
            caseDir.expand();
            caseDir.toAbsolute();
        }
    }
    else
    {
        // Nothing specified, use the current dir
        caseDir = cwd();
    }

    // The caseDir is a cleaned, absolute path

    rootPath_   = caseDir.path();
    globalCase_ = caseDir.name();
    case_       = globalCase_;  // The (processor) local case name

    // OPENFOAM API
    setEnv("FOAM_API", std::to_string(foamVersion::api), true);

    // Global case (directory) and case-name as environment variables
    setEnv("FOAM_CASE", caseDir, true);
    setEnv("FOAM_CASENAME", globalCase_, true);

    // Executable name, unless already present in the environment
    setEnv("FOAM_EXECUTABLE", executable_, false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::argList::argList
(
    int& argc,
    char**& argv,
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
:
    runControl_(),
    args_(argc),
    options_(argc),
    libs_()
{
    // Pre-scan for some options needed for initial setup:
    //   -fileHandler (takes an argument)
    //   -mpi-threads (bool option)
    //
    // Also handle -dry-run and -verbose counting
    // (it is left to the application to decide what to do with them).
    // Detect any parallel run options

    word fileHandlerName;

    if (parallelThreads_)
    {
        // Default -mpi-threads configured statically from application
        runControl_.threads(true);
    }

    for (int argi = 1; argi < argc; ++argi)
    {
        const char *optName = argv[argi];

        if (optName[0] == '-')
        {
            ++optName;  // Looks like an option, skip leading '-'
            bool emitErrorMessage = false;

            if (strcmp(optName, "dry-run") == 0)
            {
                runControl_.incrDryRun();
            }
            else if (strcmp(optName, "verbose") == 0)
            {
                runControl_.incrVerbose();
            }
            else if (strcmp(optName, "mpi-threads") == 0)
            {
                runControl_.threads(true);
            }
            else if (strcmp(optName, "fileHandler") == 0)
            {
                // Requires a parameter
                if (argi < argc-1)
                {
                    ++argi;
                    fileHandlerName = argv[argi];
                }
                else
                {
                    emitErrorMessage = true;
                }
            }
            else if (validParOptions.contains(optName))
            {
                // Contains a parallel run option
                runControl_.parRun(true);
            }

            if (emitErrorMessage)
            {
                // Missing argument: emit message but not exit or
                // FatalError since Pstream etc are not yet initialised

                Info<< nl
                    << "Error: option '-" << optName
                    << "' requires an argument" << nl << nl;

                //NO: UPstream::exit(1);  // works for serial and parallel
            }
        }
    }

    // No -fileHandler specifed, get from environment or use default
    if (fileHandlerName.empty())
    {
        fileHandlerName = Foam::getEnv("FOAM_FILEHANDLER");
        if (fileHandlerName.empty())
        {
            fileHandlerName = fileOperation::defaultFileHandler;
        }
    }

    // Parse out any additional fileHandler-specific options
    // (may alter argv list). Recover its threading requirements
    {
        auto fileOperationInit = fileOperations::fileOperationInitialise::New
        (
            fileHandlerName,
            argc,
            argv
        );

        if (fileOperationInit && fileOperationInit->needsThreading())
        {
            runControl_.threads(true);
        }
    }

    // Parallel job options detected?
    // - start parallel run (possibly filters argv as a side-effect)
    if (runControl_.parRun())
    {
        runControl_.runPar(argc, argv);
    }

    // addProfiling(argList, "argList");

    // ------------------------------------------------------------------------

    // Convert argv -> args_ and capture ( ... ) lists
    regroupArgv(argc, argv);
    commandLine_ += args_[0];

    // Set executable name immediately - useful when emitting errors.
    executable_ = fileName(args_[0]).name();

    // Check arguments and options, argv[0] was already handled
    int nArgs = 1;
    for (int argi = 1; argi < args_.size(); ++argi)
    {
        commandLine_ += ' ';
        commandLine_ += args_[argi];

        const char *optName = args_[argi].data();

        if (optName[0] == '-')
        {
            ++optName;  // Looks like an option, skip leading '-'

            if (!*optName)
            {
                Warning
                    << "Ignoring lone '-' on the command-line" << endl;
                continue;
            }

            // Option known and expects an argument?
            // - use Switch for a tri-state
            //   True  : known option, expects a parameter
            //   False : known option, no parameter
            //   bad() : unknown option

            Switch wantArg(Switch::INVALID);
            auto optIter = validOptions.cfind(optName);
            if
            (
                optIter.good()
             || (optIter = validParOptions.cfind(optName)).good()
            )
            {
                wantArg = !optIter.val().empty();
            }


            if (wantArg)
            {
                // Option expects a parameter
                // - get it or emit a FatalError.

                ++argi;
                if (argi >= args_.size())
                {
                    foamVersion::printBuildInfo(Info.stdStream(), false);

                    Info<< nl
                        <<"Error: option '-" << optName
                        << "' requires an argument" << nl << nl
                        << "See '" << executable_ << " -help' for usage" << nl
                        << " or '" << executable_
                        << " -help-full' for extended usage" << nl
                        << nl;

                    UPstream::exit(1);  // works for serial and parallel
                }

                commandLine_ += ' ';
                commandLine_ += args_[argi];

                //
                // Special handling of these options
                //

                if (strcmp(optName, "lib") == 0)
                {
                    // The '-lib' option:
                    // Add name(s) to libs for later opening
                    libs().push_back(this->getList<fileName>(argi));
                }
                else if (strcmp(optName, "debug-switch") == 0)
                {
                    // The '-debug-switch' option:
                    // change registered debug switch
                    DetailInfo << "debug-switch ";
                    debug::debugObjects()
                        .setNamedValue(args_[argi], 1, true);
                }
                else if (strcmp(optName, "info-switch") == 0)
                {
                    // The '-info-switch' option:
                    // change registered info switch
                    DetailInfo << "info-switch ";
                    debug::infoObjects()
                        .setNamedValue(args_[argi], 1, true);
                }
                else if (strcmp(optName, "opt-switch") == 0)
                {
                    // The '-opt-switch' option:
                    // change registered optimisation switch
                    DetailInfo << "opt-switch ";
                    debug::optimisationObjects()
                        .setNamedValue(args_[argi], 1, true);
                }
                else
                {
                    // Regular option (with a parameter):
                    // Duplicates handled by using the last -option specified
                    options_.set(optName, args_[argi]);
                }
            }
            else
            {
                // All other options (including unknown ones) are simply
                // registered as existing.

                options_.insert(optName, "");

                // // Special increment handling for some known flags
                // if (wantArg.good())
                // {
                //     ...
                // }
            }
        }
        else
        {
            if (nArgs != argi)
            {
                args_[nArgs] = args_[argi];
            }
            ++nArgs;
        }
    }

    args_.resize(nArgs);

    parse(checkArgs, checkOpts, initialise);
}


Foam::argList::argList
(
    const argList& args,
    const HashTable<string>& options,
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
:
    runControl_(args.runControl_),
    args_(args.args_),
    options_(options),
    libs_(),
    executable_(args.executable_),
    rootPath_(args.rootPath_),
    globalCase_(args.globalCase_),
    case_(args.case_),
    commandLine_(args.commandLine_)
{
    parse(checkArgs, checkOpts, initialise);
}


void Foam::argList::parse
(
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
{
    // Help/documentation options:
    //   -doc         Display documentation in browser
    //   -doc-source  Display source code in browser
    //   -help        Display short help and exit
    //   -help-compat Display compatibility options
    //   -help-full   Display full help and exit
    {
        bool quickExit = false;

        // Display either application or source documentation, not both
        if (options_.contains("doc"))
        {
            displayDoc(false);
            quickExit = true;
        }
        else if (options_.contains("doc-source"))
        {
            displayDoc(true);
            quickExit = true;
        }

        // Display either short or full help, not both
        if (options_.contains("help-full"))
        {
            printUsage(true);
            quickExit = true;
        }
        else if (options_.contains("help-notes"))
        {
            printNotes();
            Info<< nl;
            quickExit = true;
        }
        else if (options_.contains("help"))
        {
            printUsage(false);
            quickExit = true;
        }
        else if (options_.contains("help-man"))
        {
            printMan();
            quickExit = true;
        }

        // Allow independent display of compatibility information
        if (options_.contains("help-compat"))
        {
            printCompat();
            quickExit = true;
        }

        if (quickExit)
        {
            std::exit(0);
        }
    }

    // Print the collected error messages and exit if check fails
    if (!check(checkArgs, checkOpts))
    {
        foamVersion::printBuildInfo(Info.stdStream(), false);
        FatalError.write(Info, false);

        UPstream::exit(1);  // works for serial and parallel
    }

    if (initialise)
    {
        const string dateString = clock::date();
        const string timeString = clock::clockTime();

        // Print the banner once only for parallel runs
        if (UPstream::master() && bannerEnabled())
        {
            IOobject::writeBanner(Info, true)
                << "Build  : ";

            if (foamVersion::build.size())
            {
                Info<< foamVersion::build.c_str() << ' ';
            }

            Info<< "OPENFOAM=" << foamVersion::api;

            if (foamVersion::patched())
            {
                // Patch-level, when defined
                Info<< " patch=" << foamVersion::patch.c_str();
            }

            Info<< " version=" << foamVersion::version.c_str();

            Info<< nl
                << "Arch   : " << foamVersion::buildArch << nl
                << "Exec   : " << commandLine_.c_str() << nl
                << "Date   : " << dateString.c_str() << nl
                << "Time   : " << timeString.c_str() << nl
                << "Host   : " << Foam::hostName().c_str() << nl
                << "PID    : " << pid() << nl;
        }

        jobInfo.add("startDate", dateString);
        jobInfo.add("startTime", timeString);
        jobInfo.add("userName", userName());

        jobInfo.add("foamApi", foamVersion::api);
        jobInfo.add("foamVersion", word(foamVersion::version));

        // Add build information - only use the first word
        {
            std::string build(foamVersion::build);
            const auto space = build.find(' ');
            if (space != std::string::npos)
            {
                build.resize(space);
            }
            jobInfo.add("foamBuild", build);
        }

        jobInfo.add("code", executable_);
        jobInfo.add("argList", commandLine_);
        jobInfo.add("currentDir", cwd());
        jobInfo.add("PPID", ppid());
        jobInfo.add("PGID", pgid());

        // Load additional libraries (verbosity according to banner setting)
        libs().open(bannerEnabled());
    }


    // Set fileHandler. In increasing order of priority:
    // 1. default = uncollated
    // 2. env variable "FOAM_FILEHANDLER"
    // 3. etc/controlDict optimisationSwitches 'fileHandler'
    // 4. system/controlDict 'fileHandler' (not handled here; done in TimeIO.C)
    // 5. '-fileHandler' commmand-line option

    {
        word fileHandlerName
        (
            options_.lookup("fileHandler", Foam::getEnv("FOAM_FILEHANDLER"))
        );

        if (fileHandlerName.empty())
        {
            fileHandlerName = fileOperation::defaultFileHandler;
        }

        (void) fileOperation::fileHandler
        (
            fileOperation::New(fileHandlerName, bannerEnabled())
        );
    }


    stringList hostMachine;
    stringList hostProcs;
    const int writeHostsSwitch = Foam::debug::infoSwitch("writeHosts", 1);
    const int writeRootsSwitch = Foam::debug::infoSwitch("writeRoots", 1);

    // Collect machine/pid, and check that the build is identical
    if (runControl_.parRun())
    {
        if (UPstream::master())
        {
            hostMachine.resize(UPstream::nProcs()-1);
            hostProcs.resize(UPstream::nProcs()-1);
            string procBuild;
            label procPid;
            int proci = 0;
            for (const int subproci : UPstream::subProcs())
            {
                IPstream fromSubproc(UPstream::commsTypes::scheduled, subproci);

                fromSubproc >> procBuild >> hostMachine[proci] >> procPid;

                hostProcs[proci] = hostMachine[proci] + "." + name(procPid);
                ++proci;

                // Verify that all processors are running the same build
                if (procBuild != foamVersion::build)
                {
                    FatalErrorIn(executable())
                        << "Running build version " << foamVersion::build
                        << " but proc " << subproci << " is running "
                        << procBuild << nl
                        << exit(FatalError);
                }
            }
        }
        else
        {
            OPstream toMaster
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo()
            );
            toMaster << foamVersion::build << Foam::hostName() << Foam::pid();
        }
    }


    // Case is a single processor run unless it is running parallel
    int nProcs = 1;

    // Roots if running distributed. Only sized on the master
    fileNameList roots;

    enum distributedCodes
    {
        NON_DISTRIBUTED = 0,
        DISTRIBUTED = 1,
        DISTRIBUTED_SINGLE_ROOT = 2,
        DISTRIBUTED_MULTIPLE_ROOTS = 3
    };

    // Track which type of distributed roots etc are being used
    label distributedType
    (
        runControl_.distributed()
      ? distributedCodes::DISTRIBUTED
      : distributedCodes::NON_DISTRIBUTED
    );

    // Some cases where knowing the writeFormat can be useful...
    // label writeFormat(-1);


    // If this actually is a parallel run
    if (runControl_.parRun())
    {
        // For the master
        if (UPstream::master())
        {
            // Establish rootPath_/globalCase_/case_ for master
            setCasePaths();

            // The system/decomposeParDict (or equivalent)
            fileName source;

            if (this->readIfPresent("decomposeParDict", source))
            {
                bool adjustOpt = false;

                if (Foam::isDir(source))
                {
                    source /= "decomposeParDict";
                    adjustOpt = true;
                }

                // Case-relative if not absolute and not "./" etc
                if (!source.isAbsolute() && !source.starts_with('.'))
                {
                    source = rootPath_/globalCase_/source;
                    adjustOpt = true;
                }

                // Could also check for absolute path, but shouldn't be needed
                if (adjustOpt)
                {
                    source.clean();  // Remove unneeded ".."
                    options_.set("decomposeParDict", source);
                }
            }

            // If running distributed (different roots for different procs)
            label dictNProcs = -1;
            if (this->readListIfPresent("roots", roots))
            {
                source = "-roots";
                runControl_.distributed(true);
                distributedType = distributedCodes::DISTRIBUTED;

                if (roots.empty())
                {
                    FatalErrorInFunction
                        << "The -roots option must contain values"
                        << exit(FatalError);
                }
                if (roots.size() > 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else if (options_.contains("hostRoots"))
            {
                source = "-hostRoots";
                runControl_.distributed(true);
                distributedType = distributedCodes::DISTRIBUTED;

                ITstream is(this->lookup("hostRoots"));

                List<Tuple2<wordRe, fileName>> hostRoots(is);
                checkITstream(is, "hostRoots");

                if (hostRoots.empty())
                {
                    FatalErrorInFunction
                        << "The -hostRoots option must contain values"
                        << exit(FatalError);
                }

                // Match machine names to roots
                roots.resize(UPstream::nProcs()-1, fileName::null);
                for (const auto& hostRoot : hostRoots)
                {
                    labelList matched
                    (
                        wordRes::matching(hostRoot.first(), hostMachine)
                    );
                    for (const label matchi : matched)
                    {
                        if (!roots[matchi].empty())
                        {
                            FatalErrorInFunction
                                << "Multiple matching roots for "
                                << hostMachine[matchi] << " in "
                                << hostRoots << nl
                                << exit(FatalError);
                        }

                        roots[matchi] = hostRoot.second();
                    }
                }

                // Check
                forAll(roots, hosti)
                {
                    if (roots[hosti].empty())
                    {
                        FatalErrorInFunction
                            << "No matching roots for "
                            << hostMachine[hosti] << " in "
                            << hostRoots << nl
                            << exit(FatalError);
                    }
                }

                if (roots.size() > 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else if (checkProcessorDirectories_ && UPstream::nProcs() > 1)
            {
                // Check values from decomposeParDict

                const bool useDefault = source.empty();
                if (useDefault)
                {
                    source = rootPath_/globalCase_/"system"/"decomposeParDict";
                }

                // Disable any parallel comms happening inside the fileHandler
                // since we are on master. This can happen e.g. inside
                // the masterUncollated/collated handler. Note that we
                // also have to protect the actual dictionary parsing since
                // it might trigger file access (e.g. #include, #codeStream)
                const bool oldParRun = UPstream::parRun(false);
                // Note: non-parallel running might update
                // fileOperation::nProcs() so store & restore below
                const label nOldProcs = fileHandler().nProcs();

                autoPtr<ISstream> dictStream
                (
                    fileHandler().NewIFstream(source)
                );

                if (dictStream && dictStream->good())
                {
                    // Get numberOfSubdomains if it exists.
                    // - mandatory when running with distributed roots

                    IOobjectOption::readOption nDomainsReadOpt
                        = IOobjectOption::READ_IF_PRESENT;

                    dictionary decompDict(*dictStream);

                    if (decompDict.getOrDefault("distributed", false))
                    {
                        nDomainsReadOpt = IOobjectOption::MUST_READ;
                        runControl_.distributed(true);
                        distributedType = distributedCodes::DISTRIBUTED;

                        decompDict.readEntry("roots", roots);

                        if (roots.empty())
                        {
                            DetailInfo
                                << "WARNING: running distributed"
                                << " but did not specify roots!" << nl;
                        }
                    }

                    decompDict.readEntry
                    (
                        "numberOfSubdomains",
                        dictNProcs,
                        keyType::LITERAL,
                        nDomainsReadOpt
                    );
                }
                else
                {
                    if (useDefault)
                    {
                        // Optional if using default location
                        DetailInfo
                            << "WARNING: running without decomposeParDict "
                            << this->relativePath(source) << nl;
                    }
                    else
                    {
                        // Mandatory if specified as -decomposeParDict
                        FatalError
                            << "Cannot read decomposeParDict: "
                            << this->relativePath(source) << nl
                            << exit(FatalError);
                    }
                }

                UPstream::parRun(oldParRun);  // Restore parallel state
                const_cast<fileOperation&>(fileHandler()).nProcs(nOldProcs);

                if (UPstream::nProcs() == 1)
                {
                    Warning
                        << "Running parallel on single processor. This only"
                        << " makes sense for multi-world simulation" << endl;
                    dictNProcs = 1;
                }
            }


            // Distributed roots
            if (!roots.empty())
            {
                for (fileName& dir : roots)
                {
                    dir.expand();
                }

                // Identical root specified everywhere?
                // - use optimized single-root variant
                if (roots.size() > 1 && roots.uniform())
                {
                    roots.resize(1);
                }

                if (roots.size() == 1)
                {
                    // Single root specified, use it for all processes
                    distributedType =
                        distributedCodes::DISTRIBUTED_SINGLE_ROOT;

                    // Adjust dictNProcs for command-line '-roots' option
                    if (dictNProcs <= 0)
                    {
                        dictNProcs = UPstream::nProcs();
                    }
                }
                else if (roots.size() > 1)
                {
                    distributedType =
                        distributedCodes::DISTRIBUTED_MULTIPLE_ROOTS;

                    if (roots.size() != UPstream::nProcs()-1)
                    {
                        FatalError
                            << "Number of roots " << roots.size()
                            << " != number of sub-ranks "
                            << UPstream::nProcs()-1
                            << exit(FatalError);
                    }
                }
            }


            //
            // Check number of processors.
            //

            // nProcs     => number of actual procs
            // dictNProcs => number of procs specified in decompositionDict
            // nProcDirs  => number of processor directories
            //               (n/a when running distributed)
            //
            // - normal running : nProcs = dictNProcs = nProcDirs
            // - decomposition to more  processors : nProcs = dictNProcs
            // - decomposition to fewer processors : nProcs = nProcDirs
            if
            (
                checkProcessorDirectories_
             && UPstream::nProcs() > 1
            )
            {
                if (dictNProcs > UPstream::nProcs())
                {
                    FatalError
                        << this->relativePath(source)
                        << " specifies " << dictNProcs
                        << " processors but job was started with "
                        << UPstream::nProcs() << " ranks."
                        << exit(FatalError);
                }

                // Possibly going to fewer processors.
                // Check if all procDirs are there.
                // NOTE: Only works when not using distributed roots!
                if
                (
                    // Can only rely on directory scanning *without* distributed roots!
                    roots.empty()
                 && dictNProcs >= 1
                 && dictNProcs < UPstream::nProcs()
                )
                {
                    label nProcDirs = 0;
                    {
                        const bool oldParRun = UPstream::parRun(false);
                        // Don't cache processor directories (probably not
                        // needed since master-only
                        const int oldCacheLevel(fileOperation::cacheLevel(0));
                        // Accept any processorsXXX
                        const int oldFilter(fileOperation::nProcsFilter(0));

                        nProcDirs = fileHandler().nProcs
                        (
                            rootPath_/globalCase_,
                            "",
                            dictNProcs  // expected nProcs
                        );

                        fileOperation::nProcsFilter(oldFilter);
                        fileOperation::cacheLevel(oldCacheLevel);
                        UPstream::parRun(oldParRun);
                    }

                    if (nProcDirs < UPstream::nProcs())
                    {
                        FatalError
                            << "Number of processor directories = " << nProcDirs
                            << " is not equal to the number of ranks = "
                            << UPstream::nProcs()
                            << exit(FatalError);
                    }
                }
            }


            // Broadcast the master's argument list (unaltered)
            {
                OPBstream toProcs(UPstream::worldComm);

                toProcs
                    << args_ << options_
                    << distributedType
                    << label(runControl_.dryRun())
                    << label(runControl_.verbose());
            }
        }
        else
        {
            // Receive the broadcasted master's argument list
            label numDryRun, numVerbose;

            IPBstream fromMaster(UPstream::worldComm);

            fromMaster
                >> args_ >> options_
                >> distributedType
                >> numDryRun >> numVerbose;

            runControl_.distributed(distributedType);
            runControl_.dryRun(numDryRun);
            runControl_.verbose(numVerbose);
        }


        // Final handling of distributed roots (if any)
        if
        (
            distributedType == distributedCodes::DISTRIBUTED_SINGLE_ROOT
        )
        {
            // The same root for all sub-ranks
            // - use broadcast to transmit value

            fileName newCasePath;

            if (UPstream::master())
            {
                newCasePath = roots[0]/globalCase_;
                OPBstream::send(newCasePath);  // worldComm
            }
            else
            {
                IPBstream::recv(newCasePath);  // worldComm
                options_.set("case", newCasePath);
            }
        }
        else if
        (
            distributedType == distributedCodes::DISTRIBUTED_MULTIPLE_ROOTS
        )
        {
            // Different roots for each sub-rank
            // - use point-to-point communication to transmit values

            fileName newCasePath;

            if (UPstream::master())
            {
                for (const int subproci : UPstream::subProcs())
                {
                    newCasePath = roots[subproci-1]/globalCase_;
                    OPstream::send(newCasePath, subproci);  // worldComm
                }
            }
            else
            {
                IPstream::recv(newCasePath, UPstream::masterNo());  // worldComm
                options_.set("case", newCasePath);
            }
        }


        // Establish rootPath_/globalCase_/case_ for sub-process
        if (!UPstream::master())
        {
            setCasePaths();
        }


        nProcs = UPstream::nProcs();
        if (UPstream::nProcs() > 1)
        {
            case_ =
            (
                globalCase_
              / ("processor" + Foam::name(UPstream::myProcNo()))
            );
        }
        else
        {
            case_ = globalCase_;
        }
    }
    else
    {
        // Establish rootPath_/globalCase_/case_
        setCasePaths();
        case_ = globalCase_;   // Redundant, but extra safety?
    }


    // If needed, adjust fileHandler for distributed roots
    if (runControl_.distributed() && fileOperation::fileHandlerPtr_)
    {
        fileOperation::fileHandlerPtr_->distributed(true);

        const labelList& ranks = fileHandler().ioRanks();

        if (runControl_.parRun() && ranks.size())
        {
            // Detect processor directories both on local proc and on
            // (world) master proc. If the local proc doesn't have them
            // but the master has it will attempt to copy them.

            // Expected local directory name
            const fileName procDir
            (
                rootPath_
              / globalCase_
              / ("processor" + Foam::name(UPstream::myProcNo()))
            );

            // Try and find my local directory using the fileHandler. This
            // will check the local disk on the IO rank
            // (since running distributed)
            fileNameList pathDirs(UPstream::nProcs());
            auto& pathDir = pathDirs[UPstream::myProcNo()];
            pathDir = fileHandler().filePath(procDir, false);

            if (returnReduceOr(pathDir.empty()))
            {
                // There is at least one processor that cannot find
                // the processor directory. Look for it on the master.
                // E.g. decomposed into 4 processors, two roots:
                //      processors4_0-1/
                //      processors4_2-3/
                // So:
                // - processor0 reads the same disk as processor0
                // - processor2 needs the whole directory sent over
                // - processor3 reads the same disk as processor2
                if (UPstream::master() && bannerEnabled())
                {
                    Info<< "I/O    :"
                        << " distributed - copying missing directories"
                        << nl;
                }

                // Collect all wanted directories (or empty). Note: could
                // just collect missing ones ...

                Pstream::gatherList(pathDirs);
                fileName masterRootPath(rootPath_);
                Pstream::broadcast(masterRootPath);

                List<fileNameList> rankToDirs(UPstream::nProcs());
                if (UPstream::master())
                {
                    const bool oldParRun = UPstream::parRun(false);
                    // Note: non-parallel running might update
                    // fileOperation::nProcs() so store & restore below
                    const label nOldProcs = fileHandler().nProcs();

                    label rank = 0;
                    for (label proci = 1; proci < pathDirs.size(); ++proci)
                    {
                        if (ranks.contains(proci))
                        {
                            rank = proci;
                        }

                        if (pathDirs[proci].empty())
                        {
                            // Synthesise corresponding name on the master
                            // processor
                            const fileName procDir
                            (
                                rootPath_
                              / globalCase_
                              / ("processor" + Foam::name(proci))
                            );
                            const fileName foundDir
                            (
                                fileHandler().filePath(procDir, false)
                            );

                            if
                            (
                                !foundDir.empty()
                             && !rankToDirs[rank].contains(foundDir)
                            )
                            {
                                rankToDirs[rank].push_back(foundDir);
                            }
                        }
                    }

                    UPstream::parRun(oldParRun);
                    const_cast<fileOperation&>(fileHandler()).nProcs(nOldProcs);
                }
                Pstream::broadcast(rankToDirs);

                // Copy missing directories on all the IOranks.
                // Note: instead of passing 'writeOnProc' flag we could create
                // communicator just between master and IOrank, but that is
                // also expensive.

                forAll(rankToDirs, proci)
                {
                    // Am I the reponsible IOrank for this processor
                    const bool amIO = (UPstream::myProcNo() == proci);

                    // Construct equivalent directory on proci
                    for (const auto& srcDir : rankToDirs[proci])
                    {
                        const fileName tgtDir
                        (
                            rootPath_
                          / srcDir.relative(masterRootPath)
                        );

                        if (amIO)
                        {
                            // I am the IO rank
                            Pout<< "On rank " << proci << nl
                                << "    copying : " << srcDir << nl
                                << "    to      : " << tgtDir << endl;
                        }

                        fileHandler().broadcastCopy
                        (
                            UPstream::worldComm,
                            amIO,
                            tgtDir,
                            tgtDir
                        );
                    }
                }
            }
        }
    }


    // Keep/discard sub-process host/root information for reporting:
    if (UPstream::master() && runControl_.parRun())
    {
        if (!writeHostsSwitch)
        {
            // Clear here to ensures it doesn't show in the jobInfo
            hostProcs.clear();
        }
        if (!writeRootsSwitch)
        {
            roots.clear();
        }
    }

    if (UPstream::master() && bannerEnabled())
    {
        Info<< "Case   : " << (rootPath_/globalCase_).c_str() << nl
            << "nProcs : " << nProcs << nl;

        if (runControl_.parRun())
        {
            if (hostProcs.size())
            {
                if (writeHostsSwitch == 1)
                {
                    // Compact output (see etc/controlDict)
                    printHostsSubscription(hostProcs);
                }
                else if (writeHostsSwitch)
                {
                    // Full output of "host.pid"
                    Info<< "Hosts  :\n(" << nl;

                    // Include master in the list
                    Info<< "    " << Foam::hostName().c_str() << '.'
                        << Foam::pid() << nl;

                    // Sub-processes
                    for (const auto& str : hostProcs)
                    {
                        Info<< "    " << str.c_str() << nl;
                    }
                    Info<< ')' << nl;
                }
            }
            if (roots.size())
            {
                bool hasPrinted = false;
                if (writeRootsSwitch == 1)
                {
                    // Compact output
                    hasPrinted = printRootsSubscription(hostProcs, roots);
                }
                if (writeRootsSwitch && !hasPrinted)
                {
                    // Full output
                    Info<< "Roots  : " << roots << nl;
                }
            }
            Info<< "Pstream initialized with:" << nl
                << "    floatTransfer      : "
                << Switch::name(UPstream::floatTransfer) << nl
                << "    maxCommsSize       : "
                << UPstream::maxCommsSize << nl
                << "    nProcsSimpleSum    : "
                << UPstream::nProcsSimpleSum << nl
                << "    nonBlockingExchange: "
                << UPstream::nProcsNonblockingExchange
                << " (tuning: " << UPstream::tuning_NBX_ << ')' << nl
                << "    exchange algorithm : "
                << PstreamBuffers::algorithm << nl
                << "    commsType          : "
                << UPstream::commsTypeNames[UPstream::defaultCommsType] << nl
                << "    polling iterations : "
                << UPstream::nPollProcInterfaces << nl;

            if (UPstream::allWorlds().size() > 1)
            {
                Info<< "    worlds             : "
                    << flatOutput(UPstream::allWorlds()) << nl
                    << "    world              : " << UPstream::myWorld()
                    << nl;
            }
        }
    }

    if (initialise)
    {
        jobInfo.add("root", rootPath_);
        jobInfo.add("case", globalCase_);
        jobInfo.add("nProcs", nProcs);
        if (hostProcs.size())
        {
            jobInfo.add("hosts", hostProcs);
        }
        if (roots.size())
        {
            jobInfo.add("roots", roots);
        }
        jobInfo.write();

        // Switch on signal trapping. We have to wait until after Pstream::init
        // since this sets up its own ones.
        sigFpe::set(bannerEnabled());
        sigInt::set(bannerEnabled());
        sigQuit::set(bannerEnabled());
        sigSegv::set(bannerEnabled());

        if (UPstream::master() && bannerEnabled())
        {
            Info<< "fileModificationChecking : "
                << "Monitoring run-time modified files using "
                << IOobject::fileCheckTypesNames
                    [
                        IOobject::fileModificationChecking
                    ];
            if
            (
                IOobject::fileModificationChecking == IOobject::timeStamp
             || IOobject::fileModificationChecking == IOobject::timeStampMaster
            )
            {
                if (IOobject::maxFileModificationPolls == 1)
                {
                    Info<< " (fileModificationSkew "
                        << IOobject::fileModificationSkew
                        << ")";
                }
                else if (IOobject::maxFileModificationPolls > 1)
                {
                    Info<< " (fileModificationSkew "
                        << IOobject::fileModificationSkew
                        << ", maxFileModificationPolls "
                        << IOobject::maxFileModificationPolls
                        << ")";
                }
                else
                {
                    FatalErrorInFunction
                        << "Invalid setting for maxFileModificationPolls "
                        << IOobject::maxFileModificationPolls
                        << exit(FatalError);
                }
            }
            Info<< nl;

            Info<< "allowSystemOperations : ";
            if (dynamicCode::allowSystemOperations)
            {
                Info<< "Allowing";
            }
            else
            {
                Info<< "Disallowing";
            }
            Info<< " user-supplied system call operations" << nl
                << nl;
            IOobject::writeDivider(Info);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::argList::~argList()
{
    jobInfo.stop();     // Normal job termination

    // Delete file handler to flush any remaining IO
    (void) fileOperation::fileHandler(nullptr);
}


// * * * * * * * * * * * * * * * Capabilities  * * * * * * * * * * * * * * * //

bool Foam::argList::allowFunctionObjects() const
{
    if (validOptions.contains("withFunctionObjects"))
    {
        // '-withFunctionObjects' is available and explicitly enabled
        return options_.contains("withFunctionObjects");
    }
    else if (validOptions.contains("noFunctionObjects"))
    {
        // '-noFunctionObjects' is available and not explicitly disabled
        return !options_.contains("noFunctionObjects");
    }

    // Disallow functions if there is no way to enable/disable them
    return false;
}


bool Foam::argList::allowLibs() const
{
    return !options_.contains("no-libs");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::argList::count(const UList<word>& optionNames) const
{
    label num = 0;
    for (const word& optName : optionNames)
    {
        if (options_.contains(optName))
        {
            ++num;
        }
    }
    return num;
}


Foam::label Foam::argList::count(std::initializer_list<word> optionNames) const
{
    label num = 0;
    for (const word& optName : optionNames)
    {
        if (options_.contains(optName))
        {
            ++num;
        }
    }
    return num;
}


bool Foam::argList::setOption(const word& optName, const string& param)
{
    // Some options are always protected
    if
    (
        optName == "case"
     || optName == "parallel"
     || optName == "roots"
    )
    {
        FatalErrorInFunction
            <<"Option: '" << optName << "' is protected" << nl
            << exit(FatalError);
        return false;
    }

    const auto optIter = options_.cfind(optName);

    if (!optIter.good() || (optIter.val() != param))
    {
        options_.set(optName, param);
        return true;
    }

    return false;
}


bool Foam::argList::unsetOption(const word& optName)
{
    // Some options are always protected
    if
    (
        optName == "case"
     || optName == "parallel"
     || optName == "roots"
    )
    {
        FatalErrorInFunction
            <<"Option: '" << optName << "' is protected" << nl
            << exit(FatalError);
        return false;
    }

    // Remove the option, return true if state changed
    return options_.erase(optName);
}


void Foam::argList::displayDoc(bool source) const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    fileNameList docDirs(docDict.get<fileNameList>("doxyDocDirs"));
    fileName docExt(docDict.get<fileName>("doxySourceFileExt"));

    // For source code: change xxx_8C.html to xxx_8C_source.html
    if (source)
    {
        docExt.replace(".", "_source.");
    }

    fileName url;

    for (const fileName& dir : docDirs)
    {
        // The http protocols are last in the list
        if (dir.starts_with("http:") || dir.starts_with("https:"))
        {
            url = dir/executable_ + docExt;
            break;
        }

        fileName docFile = stringOps::expand(dir/executable_ + docExt);

        if
        (
            docFile.starts_with("file://")
          ? isFile(docFile.substr(7))   // check part after "file://"
          : isFile(docFile)
        )
        {
            url = std::move(docFile);
            break;
        }
    }

    if (url.empty())
    {
        Info<< nl
            << "No documentation found for " << executable_
            << ", but you can use -help to display the usage\n" << endl;

        return;
    }

    string docBrowser = getEnv("FOAM_DOC_BROWSER");
    if (docBrowser.empty())
    {
        docDict.readEntry("docBrowser", docBrowser);
    }

    // Can use FOAM_DOC_BROWSER='application file://%f' if required
    if (docBrowser.contains("%f"))
    {
        docBrowser.replace("%f", url);
    }
    else
    {
        docBrowser += " " + url;
    }

    // Split on whitespace to use safer version of Foam::system()

    CStringList command(stringOps::splitSpace(docBrowser));

    Info
        << "OpenFOAM " << foamVersion::api << " documentation:" << nl
        << "    " << command << nl << endl;

    Foam::system(command, true);
}


bool Foam::argList::check(bool checkArgs, bool checkOpts) const
{
    bool ok = true;

    if (UPstream::master())
    {
        const label nargs = args_.size()-1;
        if (checkArgs && nargs != validArgs.size())
        {
            FatalError
                << "Expected " << validArgs.size()
                << " arguments but found " << nargs << endl;
            ok = false;
        }

        if (checkOpts)
        {
            forAllConstIters(options_, iter)
            {
                const word& optName = iter.key();
                if
                (
                    !validOptions.contains(optName)
                 && !validParOptions.contains(optName)
                )
                {
                    FatalError
                        << "Invalid option: -" << optName << endl;
                    ok = false;
                }
            }
        }

        if (!ok)
        {
            FatalError
                << nl
                << "See '" << executable_ << " -help' for usage" << nl
                << " or '" << executable_
                << " -help-full' for extended usage" << nl
                << nl;
        }
    }

    return ok;
}


bool Foam::argList::checkRootCase() const
{
    if (!fileHandler().isDir(rootPath()))
    {
        FatalError
            << executable_
            << ": cannot open root directory " << rootPath()
            << endl;

        return false;
    }

    const fileName pathDir(fileHandler().filePath(path(), false));

    if (checkProcessorDirectories_ && pathDir.empty() && UPstream::master())
    {
        // Allow non-existent processor directories on sub-processes,
        // to be created later (e.g. redistributePar)
        FatalError
            << executable_
            << ": cannot open case directory " << path()
            << endl;

        return false;
    }

    return true;
}


// ************************************************************************* //
