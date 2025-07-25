/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Test-dictionary3

Description
    Test expressions and re-expansions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "dictionary.H"
#include "vector.H"
#include "SpanStream.H"
#include "formattingEntry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

    {
        ICharStream is
        (
            "value   10;"
            "scalar1 $value;"
            "scalar2 -$value;"

            // Use #eval expansion entirely
            "vector1 ${{vector($value, -$value, $value)}};"
            "vector2 ($value -$value $value);"
        );

        dictionary dict(is);

        // Add some more entries
        {
            label idx = 0;
            dictionary subdict;

            subdict.add("key", 100);

            // subdict.comment("this would be cool!");

            subdict.add
            (
                new formattingEntry(++idx, "// comment - without newline.")
            );

            subdict.add
            (
                // NB newline must be part of the content!
                new formattingEntry(++idx, "// some comment - with newline?\n")
            );

            subdict.add
            (
                // NB newline must be part of the content!
                new formattingEntry(++idx, "/* other comment */\n")
            );

            // Other - invisible
            subdict.add(new formattingEntry(++idx, token(123), false));

            // Other - visible (probably not what anyone wants!)
            subdict.add(new formattingEntry(++idx, token(456)));

            subdict.add("val", 42);

            Info<< "subdict keys:" << flatOutput(subdict.toc()) << nl;
            dict.add("subdict", std::move(subdict));
        }

        Info<< "input dictionary:" << nl;
        IOobject::writeDivider(Info);
        dict.write(Info, false);
        IOobject::writeDivider(Info);
        Info << nl;

        Info<< "value: " << dict.get<scalar>("value") << nl;

        Info<< "scalar1: " << dict.get<scalar>("scalar1") << nl;
        Info<< "scalar2: " << dict.get<scalar>("scalar2") << nl;

        Info<< "vector1: " << dict.get<vector>("vector1") << nl;
        Info<< "vector2: " << dict.get<vector>("vector2") << nl;
    }

    return 0;
}


// ************************************************************************* //
