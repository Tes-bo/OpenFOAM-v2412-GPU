/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "MinMax.H"
#include "Switch.H"

using namespace Foam;


void predicateTests_label(const word& optName, const argList& args)
{
    Info<< "predicate tests for " << optName << nl;

    const bool oldThrowingError = FatalError.throwing(true);
    try
    {
        label val;
        val = args.getCheck<label>(optName, labelMinMax::ge(0));
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError "
            << err << nl << endl;
    }

    FatalError.throwing(oldThrowingError);
}


void predicateTests_scalar(const word& optName, const argList& args)
{
    Info<< "predicate tests for " << optName << nl;

    const bool oldThrowingError = FatalError.throwing(true);
    try
    {
        scalar val;
        val = args.getCheck<scalar>(optName, scalarMinMax::ge(0));
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError "
            << err << nl << endl;
    }

    FatalError.throwing(oldThrowingError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories(); // parallel OK, but without checks

    // argList::noFunctionObjects();
    argList::addOption("label",  "value", "Test parsing of label");
    argList::addOption("scalar", "value", "Test parsing of scalar");
    argList::addOption("string", "value", "Test string lookup");
    argList::addOption("relative", "PATH", "Test relativePath");
    argList::addBoolOption
    (
        "predicates",
        "Apply some predicate tests (for label and scalar)"
    );

    // These are actually lies (never had -parseLabel, -parseScalar etc),
    // but good for testing...

    // Emits warning about it being old
    argList::addOptionCompat("label", {"parseLabel", 1612});

    // Specifying version=0 to use alias without any warnings
    argList::addOptionCompat("scalar", {"parseScalar", 0});

    // Fake a future option...
    argList::addOptionCompat("label", {"parse-label", 2112});

    // Ignore an old bool option
    argList::ignoreOptionCompat({"xml", 1700}, false);

    // Ignore an old option with arg. Specified version=0 to suppress warnings
    argList::ignoreOptionCompat({"format", 0}, true);

    // Ignore a future option? Fairly pointless
    argList::ignoreOptionCompat({"ascii", 2112}, false);

    argList::addArgument("label");
    argList::addArgument("...");
    argList::addArgument("label");
    argList::noMandatoryArgs();

    argList::addDryRunOption("Just for testing");
    argList::addVerboseOption("Increase verbosity");

    // Check -verbose before initialisation
    UPstream::debug = argList::verbose(argc, argv);

    #include "setRootCase.H"

    Pout<< "command-line ("
        << args.options().size() << " options, "
        << args.args().size() << " args)" << nl
        << "    " << args.commandLine().c_str() << nl << nl;

    Pout<< "rootPath:   " << args.rootPath() << nl
        << "globalCase: " << args.globalCaseName() << nl
        << "globalPath: " << args.globalPath() << nl
        << nl;

    Pout<< "dry-run: " << args.dryRun()
        << " verbose: " << args.verbose() << nl;


    if (args.found("relative"))
    {
        Pout<< "input path: " << args["relative"] << nl
            << "relative  : " << args.relativePath(args["relative"], true) << nl
            << nl;
    }

    Info<< "have: "
        << args.count({"label", "scalar"}) << " options" << nl;

    label ival;
    scalar sval;

    Info<< nl;

    Info<< "-label = " << flush;
    if (args.readIfPresent("label", ival))
    {
        Info<< ival << nl;

        if (args.found("predicates"))
        {
            predicateTests_label("label", args);
        }
    }
    else
    {
        Info<< "not specified" << nl;
    }

    Info<< "-scalar = " << flush;
    if (args.readIfPresent("scalar", sval))
    {
        Info<< sval << nl;

        if (args.found("predicates"))
        {
            predicateTests_scalar("scalar", args);
        }
    }
    else
    {
        Info<< "not specified" << nl;
    }


    // Using direct reading
    Info<< nl;
    if (args.found("label"))
    {
        Info<< "-label = " << args.get<label>("label")
            << " or " << args.opt<label>("label")
            #ifdef Foam_argList_1712
            << " or " << args.optionRead<label>("label")  // old-compat
            #endif
            << " or " << readLabel(args["label"])         // with function
            << nl;
    }

    if (args.found("scalar"))
    {
        Info<< "-scalar = " << args.get<scalar>("scalar")
            << " or " << args.opt<scalar>("label")
            #ifdef Foam_argList_1712
            << " or " << args.optionRead<scalar>("scalar") // old-compat
            #endif
            << " or " << readScalar(args["scalar"])        // with function
            << nl;
    }

    if (args.found("string"))
    {
        Info<< "-string = " << args.get("string")
            << " or " << args.opt("string")
            #ifdef Foam_argList_1712
            << " or " << args.optionRead<scalar>("string")  // old-compat
            #endif
            << nl;
    }


    // Arg reading
    Info<< nl;
    for (label argi=1; argi < args.size(); ++argi)
    {
        Info<< "arg[" << argi << "] = " << args.get<string>(argi)
            #ifdef Foam_argList_1712
            << " or " << args.read<label>(argi)     // old-compat
            << " or " << args.argRead<label>(argi)  // old-compat
            #endif
            << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
