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

Application
    Test-predicates

Description
    Simple tests using predicates

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "labelList.H"
#include "wordList.H"
#include "predicates.H"
#include "FlatOutput.H"
#include "regExp.H"

using namespace Foam;


template<class ListType, class UnaryPredicate>
label printMatching(const ListType& list, const UnaryPredicate& pred)
{
    label count = 0;

    Info<< "(";

    for (const auto& val : list)
    {
        if (pred(val))
        {
            if (count) Info<< ' ';
            Info<< val;
            ++count;
        }
    }

    Info<< ") => " << count << nl;

    return count;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    wordList words
    {
        "abc",
        "def",
        "hij",
        "abc_",
        "def_",
        "hij_",
    };

    labelList values(identity(40, -10));

    Info<< "words:  " << flatOutput(words) << endl;
    Info<< "values: " << flatOutput(values)  << endl;

    regExp matcher(".*_.*");

    Info<< "With '_': ";
    printMatching(words, matcher);

    Info<< "All: ";
    printMatching(words, predicates::always());

    Info<< "None: ";
    printMatching(words, predicates::never());

    Info<< "Neg values: ";
    printMatching(values, [](const label v) { return v < 0; });

    Info<< "Even values: ";
    printMatching(values, [](const label v) { return !(v % 2); });

    Info<< "All: ";
    printMatching(values, predicates::always());

    Info<< "None: ";
    printMatching(values, predicates::never());

    Info<< "zero:";
    Info<< static_cast<label>(Foam::zero{}(values)) << nl;

    return 0;
}

// ************************************************************************* //
