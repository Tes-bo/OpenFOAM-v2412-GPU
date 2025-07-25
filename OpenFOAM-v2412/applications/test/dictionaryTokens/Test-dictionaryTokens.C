/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    printDictionary

Description

    Test dictionaryTokens

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"

#include "dictionaryTokens.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addBoolOption("info", "report token info");
    argList::addBoolOption("value", "report token value");

    argList::addArgument("dict .. dictN");
    argList args(argc, argv, false, true);

    const bool optInfo  = args.found("info");
    const bool optValue = args.found("value");

    for (label argi=1; argi < args.size(); ++argi)
    {
        IFstream is(args.get<fileName>(argi));

        dictionary dict(is);

        dictionaryTokens dictTokens(dict);

        while (dictTokens.good())
        {
            if (optInfo)
            {
                // Token info
                Info<< (*dictTokens).info() << nl;
            }
            else if (optValue)
            {
                // Token value
                Info<< *dictTokens << nl;
            }
            else
            {
                // Token type
                Info<< (*dictTokens).name() << nl;
            }
            ++dictTokens;
        }

        Info<< nl;
    }

    return 0;
}


// ************************************************************************* //
