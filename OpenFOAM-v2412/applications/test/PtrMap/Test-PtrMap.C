/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include <iostream>
#include "PtrMap.H"

using namespace Foam;

template<class T>
void printTable(const PtrMap<T>& table)
{
    Info<< table.size() << nl << "(" << nl;

    forAllConstIters(table, iter)
    {
        const T* ptr = iter.val();
        Info<< iter.key() << " = ";
        if (ptr)
        {
            Info<< *ptr << " (" << name(ptr) << ")";
        }
        else
        {
            Info<< "nullptr";
        }
        Info<< nl;
    }

    Info<< ")" << endl;

    // Values only, with for-range
    Info<< "values (";
    for (auto val : table)
    {
        Info<< ' ';
        if (val)
        {
            Info<< *val;
        }
        else
        {
            Info<< "nullptr";
        }
    }
    Info<< " )" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    PtrMap<double> myTable;
    myTable.set(1, new double(42.1));
    myTable.set(2, nullptr);
    myTable.set(3, new double(3.14159));
    myTable.set(4, new double(2.718282));
    myTable.set(4, new double(1.414214));

    // Info<< myTable << endl;
    printTable(myTable);

    PtrMap<double> copy(myTable);

    // Info<< copy << endl;
    printTable(copy);
    Info<< copy << endl;

    Info<<"\nerase some existing and non-existing entries" << nl;

    auto iter = myTable.find(3);
    myTable.erase(iter);

    iter = myTable.find(1000);  // unknown key
    myTable.erase(iter);

    myTable.erase(1);
    iter = myTable.find(100000);  // unknown key

    printTable(myTable);

    PtrMap<double> moved(std::move(copy));

    Info<< nl << "test movable" << nl;
    Info<<"input:" << nl;
    printTable(copy);

    Info<<"output:" << nl;
    printTable(moved);

    PtrMap<double> other;

    Info<<"move assign" << nl;

    other = std::move(moved);
    printTable(other);

    Info<<"old" << nl;
    printTable(moved);

    return 0;
}

// ************************************************************************* //
