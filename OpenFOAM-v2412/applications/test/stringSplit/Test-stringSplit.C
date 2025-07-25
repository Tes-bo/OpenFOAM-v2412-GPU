/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Test-stringSplit

Description
    Test string splitting

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "stringOps.H"
#include "Switch.H"

using namespace Foam;

// Simple utility
template<class StringType>
void printSubStrings
(
    const StringType& str,
    const SubStrings<StringType>& split
)
{
    Info<< "string {" << str.size() << " chars} = " << str << nl
        << split.size() << " elements {" << split.length() << " chars}"
        << nl;

    unsigned i = 0;
    for (const auto s : split)
    {
        Info<< "[" << i++ << "] {" << s.length() << " chars} = "
            << s.str() << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noMandatoryArgs();
    argList::addArgument("string .. stringN");

    argList::addOption
    (
        "any",
        "delimChars",
        "test split on any delimiter characters"
    );
    argList::addOption
    (
        "sub",
        "string",
        "test split on substring"
    );
    argList::addOption
    (
        "char",
        "delim",
        "test split on specified delimiter character"
    );
    argList::addOption
    (
        "fixed",
        "int",
        "test split on fixed width"
    );
    argList::addOption
    (
        "begin",
        "int",
        "begin offset for splits"
    );
    argList::addBoolOption
    (
        "slash",
        "test split on slash (default)"
    );
    argList::addBoolOption
    (
        "space",
        "test split on space"
    );
    argList::addBoolOption
    (
        "empty",
        "preserve empty strings in split"
    );

    argList args(argc, argv);

    if (args.size() <= 1 && args.options().empty())
    {
        args.printUsage();
    }

    const label beginOffset = args.getOrDefault<label>("begin", 0);

    const bool keepEmpty = args.found("empty");

    Info<< "begin offset: " << beginOffset << nl;
    Info<< "keep empty  : " << Switch::name(keepEmpty) << nl;

    const label nopts =
        args.count({"any", "slash", "space", "sub", "fixed", "char"});


    if (args.found("any"))
    {
        const std::string& str = args["any"];
        Info<< "split on any chars" << nl
            << "=" << str << nl
            << "~~~~~~~~~~~~~~~" << nl;

        for (label argi=1; argi < args.size(); ++argi)
        {
            auto split = stringOps::splitAny(args[argi], str, beginOffset);
            printSubStrings(args[argi], split);
        }

        if (nopts == 1)
        {
            return 0;
        }
    }

    if (args.found("sub"))
    {
        const std::string& str = args["sub"];
        Info<< "split on substring" << nl
            << "=" << str << nl
            << "~~~~~~~~~~~~~~~" << nl;

        for (label argi=1; argi < args.size(); ++argi)
        {
            auto split = stringOps::split(args[argi], str, beginOffset);
            printSubStrings(args[argi], split);
        }

        if (nopts == 1)
        {
            return 0;
        }
    }

    if (args.found("space"))
    {
        Info<< "split on space" << nl
            << "~~~~~~~~~~~~~~" << nl;

        for (label argi=1; argi < args.size(); ++argi)
        {
            auto split = stringOps::splitSpace(args[argi], beginOffset);
            printSubStrings(args[argi], split);

            Info<< "pop_front(2)" << nl;
            split.pop_front(2);
            printSubStrings(args[argi], split);
        }

        if (nopts == 1)
        {
            return 0;
        }
    }

    if (args.found("char"))
    {
        const char delim = args["char"][0];

        Info<< "split on char=" << delim << nl
            << "~~~~~~~~~~~~~~" << nl;

        for (label argi=1; argi < args.size(); ++argi)
        {
            auto split =
                stringOps::split(args[argi], delim, beginOffset, keepEmpty);
            printSubStrings(args[argi], split);
        }

        if (nopts == 1)
        {
            return 0;
        }
    }

    if (args.found("fixed"))
    {
        const label width = args.get<label>("fixed");

        Info<< "split on fixed width = " << width << nl
            << "~~~~~~~~~~~~~~" << nl;

        for (label argi=1; argi < args.size(); ++argi)
        {
            auto split = stringOps::splitFixed(args[argi], width, beginOffset);
            printSubStrings(args[argi], split);
        }

        if (nopts == 1)
        {
            return 0;
        }
    }

    // Default
    if (!nopts || args.found("slash"))
    {
        const char delim = '/';

        Info<< "split on slash" << nl
            << "~~~~~~~~~~~~~~" << nl;

        for (label argi=1; argi < args.size(); ++argi)
        {
            auto split =
                stringOps::split(args[argi], delim, beginOffset, keepEmpty);
            printSubStrings(args[argi], split);
        }
    }

    return 0;
}


// ************************************************************************* //
