/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Test scalar ranges
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "scalarRanges.H"

using namespace Foam;


template<class HTable>
Ostream& printHash(const HTable& h)
{
    const auto keys = h.sortedToc();

    Info<<"keys: " << keys.size() << "(";
    for (const auto& k : keys)
    {
        Info<<" " << k;
    }
    Info<< " )" << nl;

    Info<<"vals: " << keys.size() << "(";
    for (const auto& k : keys)
    {
        Info<<" " << h[k];
    }
    Info<< " )" << nl;

    return Info;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addArgument("range,range range,range");
    argList::addNote
    (
        "Test parsing of comma-separated ranges\n\n"
        "Eg, ':10,20,40:70,1000:'"
    );

    argList args(argc, argv, false, true);

    HashTable<scalar> testHash;
    for (label i = 0; i < 26; ++i)
    {
        word key("X");
        key[0] = 'a' + i;
        testHash.insert(key, i);

        key[0] = 'A' + i;
        testHash.insert(key, 100 + i);
    }

    printHash(testHash)
        << nl;

    for (int argi = 1; argi < argc; ++argi)
    {
        Info<<"tokenize " << args[argi] << nl;

        scalarRanges ranges(args[argi]);
        Info<<"=> " << ranges << endl;

        HashTable<scalar> subHash(testHash);
        subHash.filterValues(ranges);

        Info<<"filtered:" << nl;

        printHash(subHash)
            << nl;

        subHash.filterValues(scalarRange::le(8));
        Info<<"filtered <= 8:" << nl;

        printHash(subHash)
            << nl;
    }

    return 0;
}

// ************************************************************************* //
