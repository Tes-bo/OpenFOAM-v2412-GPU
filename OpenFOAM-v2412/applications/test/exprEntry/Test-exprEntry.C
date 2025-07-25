/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Test-exprEntry

Description
    Read in the given dictionaries and attempt to use exprEntry expansion
    on any strings.

Note
   Since this is only for testing purposes, only handles simple dictionary
   entries without attempting to descend into sub-dicts.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"
#include "stringOps.H"
#include "exprString.H"

using namespace Foam;

bool hasStrings(const primitiveEntry& e)
{
    for (const token& tok : e.stream())
    {
        if (tok.isString())
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addArgument("dict .. dictN");
    argList args(argc, argv, false, true);

    if (args.size() <= 1)
    {
        Info<< "Must supply a dictionary name!" << nl;
    }

    for (label argi=1; argi < args.size(); ++argi)
    {
        IOobject::writeDivider(Info);

        IFstream is(args.get<fileName>(argi));

        const dictionary dict(is);

        Info<< "Input dictionary:" << dict << nl
            << "With any expansions" << nl << endl;

        for (const entry& dEntry : dict)
        {
            const auto* eptr = isA<primitiveEntry>(dEntry);
            if (!eptr || !hasStrings(*eptr))
            {
                continue;
            }

            const primitiveEntry& e = *eptr;
            Info<< e << endl;

            for (const token& t : e.stream())
            {
                if (t.isString())
                {
                    string str(t.stringToken());

                    const bool oldThrowingError = FatalError.throwing(true);
                    const bool oldThrowingIOErr = FatalIOError.throwing(true);

                    try
                    {
                        // Can get an error if we have things like
                        // ${{ ... $[...] }}
                        Info<< "str : " << stringOps::expand(str, dict) << nl;
                    }
                    catch (const Foam::error& err)
                    {
                        Info<< err.message().c_str() << nl;
                    }

                    try
                    {
                        // Should not trigger any errors
                        auto expr = expressions::exprString::toExpr(str, dict);

                        Info<< "expr: " << expr << nl;
                    }
                    catch (const Foam::error& err)
                    {
                        Info<< err.message().c_str() << nl;
                    }

                    FatalError.throwing(oldThrowingError);
                    FatalIOError.throwing(oldThrowingIOErr);
                    Info<< nl;
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
