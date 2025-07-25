/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

Description
    Test the tokenizing of various things
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "SpanStream.H"
#include "cpuTime.H"
#include "DynamicList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addArgument("string .. stringN");
    argList::addOption("file", "name");
    argList::addOption("repeat", "count");
    argList::addVerboseOption("report each repeat");

    argList args(argc, argv, false, true);

    const label repeat = args.getOrDefault<label>("repeat", 1);

    cpuTime timer;
    for (label count = 0; count < repeat; ++count)
    {
        const bool verbose = (args.verbose() || count == 0);

        for (label argI=1; argI < args.size(); ++argI)
        {
            const string& rawArg = args[argI];
            if (verbose)
            {
                Info<< "input string: " << rawArg << nl;
            }

            ISpanStream is(rawArg);

            DynamicList<token> tokens;

            while (is.good())
            {
                token tok(is);
                // char ch;
                // is.get(ch);
                // is.putback(ch);
                int lookahead = is.peek();

                if (verbose)
                {
                    Info<< "token: " << tok.info()
                        << "  lookahead: '" << char(lookahead) << "'"
                        << endl;
                }

                if (tok.good())
                {
                    tokens.push_back(std::move(tok));
                    if (verbose)
                    {
                        Info<< "after append: " << tok.info() << endl;
                    }
                }
            }

            if (verbose)
            {
                Info<< nl;
                IOobject::writeDivider(Info);

                Info<< "tokenList:" << tokens << endl;
            }
        }
    }

    Info<< "tokenized args " << repeat << " times in "
        << timer.cpuTimeIncrement() << " s\n\n";

    fileName inputFile;
    if (args.readIfPresent("file", inputFile))
    {
        IFstream is(inputFile);

        for (label count = 0; count < repeat; ++count)
        {
            const bool verbose = (args.verbose() || count == 0);
            label nTokens = 0;

            if (count)
            {
                is.rewind();
            }

            Info<< nl
                << "tokenizing file (pass #" << (count+1) << ") "
                << inputFile << nl
                << "state: " << is.info() << endl;

            while (is.good())
            {
                token tok(is);
                if (verbose)
                {
                    Info<< "token: " << tok.info() << endl;
                }
                ++nTokens;
            }

            if (verbose)
            {
                Info<< nl;
                IOobject::writeDivider(Info);
            }

            Info<<"pass #" << (count+1)
                << " extracted " << nTokens << " tokens" << endl;
        }

        Info<< "tokenized file " << repeat << " times in "
            << timer.cpuTimeIncrement() << " s\n\n";
    }

    return 0;
}

// ************************************************************************* //
