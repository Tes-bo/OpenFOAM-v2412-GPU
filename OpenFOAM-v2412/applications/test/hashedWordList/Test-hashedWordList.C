/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "IOstreams.H"
#include "ITstream.H"
#include "FlatOutput.H"
#include "hashedWordList.H"

using namespace Foam;

Ostream& printInfo(const hashedWordList& list, bool withAddr=false)
{
    Info<< flatOutput(list) << nl << list.lookup() << nl;
    if (withAddr)
    {
        Info<< "addr=" << name(list.cdata()) << nl;
    }

    return Info;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Test hashedWordList" << nl;

    hashedWordList list1
    {
        "this",
        "is",
        "a",
        "list",
        "of",
        "words",
    };

    Info<<nl << "From initializer_list" << nl;
    printInfo(list1, true);

    list1.sort();

    Info<<nl << "sorted" << nl;
    printInfo(list1, true);

    // Copy construct
    hashedWordList list2(list1);

    Info<<nl << "Copy construct" << nl;
    Info<<"list1: ";
    printInfo(list1, true);
    Info<<"list2: ";
    printInfo(list2, true);

    // Move construct
    hashedWordList list3(std::move(list1));

    Info<<nl << "Move construct" << nl;
    Info<<"list1: ";
    printInfo(list1, true);
    Info<<"list3: ";
    printInfo(list3, true);

    // Move assign
    list1 = std::move(list3);

    Info<<nl << "Move assign" << nl;
    Info<<"list1: ";
    printInfo(list1, true);
    Info<<"list3: ";
    printInfo(list3, true);

    list1.swap(list3);
    Info<<nl << "Swap" << nl;
    Info<<"list1: ";
    printInfo(list1, true);
    Info<<"list3: ";
    printInfo(list3, true);

    wordList wlist1
    {
        "plain", "list", "with", "some", "with", "list", "duplicates"
    };


    // Copy construct unique
    hashedWordList list4(wlist1, true);

    Info<<nl << "Copy construct unique" << nl;
    Info<<"words: " << flatOutput(wlist1) << nl;
    Info<<"list4: ";
    printInfo(list4, false);

    // Move construct unique
    hashedWordList list5(std::move(wlist1), true);

    Info<<nl << "Move construct unique" << nl;
    Info<<"words: " << flatOutput(wlist1) << nl;
    Info<<"list5: ";
    printInfo(list5, false);

    // Move back. Leaves lookup() with some rubbish, but clean that later
    wlist1 = std::move(list5);

    Info<<nl << "Move to wordList" << nl;
    Info<<"words: " << flatOutput(wlist1) << nl;
    Info<<"list5: ";
    printInfo(list5, false);

    // Move back. Leaves lookup() with some rubbish, but clean that later
    list5 = std::move(wlist1);

    Info<<nl << "Moved from wordList" << nl;
    Info<<"words: " << flatOutput(wlist1) << nl;
    Info<<"list5: ";
    printInfo(list5, false);


    // Test access:

    Info<<nl << "Access" << nl;
    Info<<"list: " << flatOutput(list5) << nl;

    for (const auto str : { "some", "list", "of", "words" })
    {
        Info<<"index of " << str << " = " << list5[str] << nl;
    }

    // Stream construct
    {
        ITstream input
        (
            "(plain list with some with list duplicates)"
        );

        hashedWordList list6(input);

        Info<<nl << "Construct from stream" << nl;
        Info<<"list: " << flatOutput(list6) << nl;

        input.rewind();

        input >> list4;
        Info<<nl << "re-read from stream" << nl;
        Info<<"list: " << flatOutput(list4) << nl;
    }


    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //
