/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Test-edges

Description
    Simple tests for edges

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "edgeList.H"
#include "edgeHashes.H"

using namespace Foam;

void printInfo(const edge& e)
{
    Info<< "edge: " << e << " count:" << e.count() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    edge e1;
    printInfo(e1);
    Info<<"has '2'? " << e1.found(2) << endl;

    edge e2(1, 2);
    printInfo(e2);
    Info<<"has '2'? " << e2.found(2) << endl;

    edge e3{2, 3};
    printInfo(e3);
    Info<<"has '2'? " << e3.found(2) << endl;

    edge e4(4, 4);
    printInfo(e4);
    Info<<"has '2'? " << e4.found(2) << endl;

    Info<<"collapse? -> " << e4.collapse() << endl;
    printInfo(e4);

    Info<< e3 << " connects " << e2 << " => " << e2.connected(e3) << endl;

    labelPair labels(e3);

    Info<< "as labelPair: " << labels << endl;

    edge e5;
    // Good: this fails (explicit constructor):   printInfo(labels);
    // Good: this also fails (no assignment operator): e5 = labels;

    // OK: explicit
    edge e6(labels);

    Info<< nl << "hash-like functionality" << nl;

    e4.clear();

    printInfo(e4);
    for (label i : {2, -1, 2, 1, 4, 1, 2, 3})
    {
        bool ok = e4.insert(i);
        Info<< "insert(" << i << ") = " << ok << " resulting ";
        printInfo(e4);
    }

    e4.clear();
    Info<< "insert from list\n";
    labelHashSet newIndices({2, -1, 2, 1, 4, 1, 2, 3});
    e4.insert(newIndices.toc());
    printInfo(e4);

    e4.clear();
    Info<< "insert from list\n";
    e4.insert({0, 5, 2, -1, 2, 1, 4, 1, 2, 3});
    printInfo(e4);

    FixedList<label, 8> otherIndices{12, 2, -1, 1, 4, 1, 2, 3};
    e4.clear();
    Info<< "insert from list: " << otherIndices << nl;
    e4.insert(otherIndices);
    printInfo(e4);

    e4.a() = e4.b();
    Info<< "erase from list: " << otherIndices << nl;
    Info<< "removed " << e4.erase(otherIndices) << " values" << nl;
    printInfo(e4);

    for (label i : {-1, 0, 1, 3})
    {
        bool ok = e4.erase(i);
        Info<< "erase(" << i << ") = " << ok << " resulting ";
        printInfo(e4);
    }

    for (label i : {-1, 0, 1, 3})
    {
        bool ok = e4.insert(i);
        Info<< "insert(" << i << ") = " << ok << " resulting ";
        printInfo(e4);
    }
    e4.flip();
    Info<< "flipped ";
    printInfo(e4);

    for (label i : {-1, 0, 1, 3})
    {
        bool ok = e4.erase(i);
        Info<< "erase(" << i << ") = " << ok << " resulting ";
        printInfo(e4);
    }

    e4.sort();
    Info<< "sorted ";
    printInfo(e4);

    return 0;
}


// ************************************************************************* //
