/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Test-splitFunctionArgs

Description
    Test splitting of function name args

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"
#include "stringOps.H"
#include "Tuple2.H"

using namespace Foam;

// Split out function name and any arguments
// - use as per functionObjectList
void testFunctionNameAndArgsSplit(const std::string& line)
{
    word funcName;
    wordRes args;
    List<Tuple2<word, string>> namedArgs;

    const auto lbracket = line.find('(');
    if (lbracket == std::string::npos)
    {
        funcName = word::validate(line);
        // No args
    }
    else
    {
        funcName = word::validate(line.substr(0, lbracket));
        std::string params;

        const auto rbracket = line.rfind(')');
        if (rbracket != std::string::npos && lbracket < rbracket)
        {
            params = line.substr(lbracket+1, (rbracket - lbracket - 1));
        }
        else
        {
            params = line.substr(lbracket+1);
        }

        Info<<"parsing: " << params << nl;

        stringOps::splitFunctionArgs(params, args, namedArgs);
    }

    Info<< nl
        << line << nl
        << "function: <" << funcName << '>' << nl
        << "    args: " << args << nl
        << "   named: " << namedArgs << nl;
}


// Split out any arguments
void testArgsSplit(const std::string& line)
{
    wordRes args;
    List<Tuple2<word, string>> namedArgs;
    stringOps::splitFunctionArgs(line, args, namedArgs);

    Info<< nl
        << line << nl
        << "    args: " << args << nl
        << "   named: " << namedArgs << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addArgument("file1 .. fileN");
    argList args(argc, argv, false, true);

    if (args.size() <= 1)
    {
        InfoErr<< "Provide a file or files to test" << nl;
    }
    else
    {
        for (label argi=1; argi < args.size(); ++argi)
        {
            IOobject::writeDivider(Info);

            const auto inputFile = args.get<fileName>(argi);
            IFstream is(inputFile);

            string line;
            while (is.getLine(line))
            {
                if (line.empty() || line[0] == '#')
                {
                    continue;
                }

                if (line.starts_with("function:"))
                {
                    auto trim = line.find(':');
                    ++trim;
                    while (isspace(line[trim]))
                    {
                        ++trim;
                    }

                    line.erase(0, trim);
                    testFunctionNameAndArgsSplit(line);
                }
                else
                {
                    testArgsSplit(line);
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
