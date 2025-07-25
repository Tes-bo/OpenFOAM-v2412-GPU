/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "argList.H"
#include "Fstream.H"
#include "ListOps.H"
#include "IndirectList.H"
#include "labelIndList.H"
#include "SortList.H"
#include "Random.H"
#include <functional>

using namespace Foam;

template<class ListType>
void printInfo(const ListType& lst)
{
    Info<< "full: " << flatOutput(lst.values()) << nl
        << "addr: " << flatOutput(lst.addressing()) << nl
        << "list: " << flatOutput(lst) << nl
        << endl;

    Info<<"for-range :";
    for (const auto& val : lst)
    {
        Info<< " " << val;
    }
    Info<< nl;
}


template<class T, class ListType>
void testFind(const T& val, const ListType& lst)
{
    Info<< nl
        << "Search for "<< val << " in " << flatOutput(lst) << nl
        <<" found() = " << lst.found(val)
        <<" find() = " << lst.find(val)
        <<" rfind() = " << lst.rfind(val)
        <<" find(2) = " << lst.find(val, 2)
        <<" rfind(2) = " << lst.rfind(val, 2) << nl
        << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    argList::addOption
    (
        "binary",
        "file",
        "write lists in binary to specified file"
    );

    argList args(argc, argv);

    List<label> completeList(20);

    forAll(completeList, i)
    {
        completeList[i] = 10*i;
    }

    Info<< "raw : " << flatOutput(completeList) << nl << endl;

    List<label> addresses{1, 0, 3, 7, 4, 8, 5, 1, 0, 3, 7, 4, 8, 5, };

    labelIndList idl1(completeList, addresses);

    printInfo(idl1);

    Info<< "list() = ";
    idl1.list().writeList(Info, 0) << endl;

    for (const label val : { 10, 30, 40, 50, 90, 80, 120 } )
    {
        testFind(val, idl1);
    }

    {
        auto iter = idl1.cbegin();
        const auto endIter = idl1.cend();

        while (iter != endIter)
        {
            // No post-fix increment:
            Info<< *iter << nl;
            ++iter;
        }
    }

    inplaceReverseList(addresses);

    idl1.addressing() = std::move(addresses);

    printInfo(idl1);

    // Test copying
    labelUIndList uidl1(idl1);
    labelIndList idl2(uidl1);
    labelIndList idl3(idl2);

    printInfo(uidl1);

    idl1.addressing().clear();
    // idl2.addressing().clear();

    Info<<"after reset addressing:" << nl << endl;

    printInfo(uidl1);
    printInfo(idl1);
    printInfo(idl2);
    printInfo(idl3);

    fileName binaryOutput;
    if (args.readIfPresent("binary", binaryOutput))
    {
        Info<<"Writing output to " << binaryOutput << endl;

        OFstream os(binaryOutput, IOstreamOption::BINARY);

        os.writeEntry("idl1", idl1);
        os.writeEntry("idl2", idl2);
        os.writeEntry("idl3", idl3);
    }

    if (UPstream::parRun())
    {
        if (UPstream::master())
        {
            Pout<< "full: " << flatOutput(idl3.values()) << nl
                << "send: " << flatOutput(idl3) << endl;

            for (const int proci : UPstream::subProcs())
            {
                OPstream::send(idl3, proci);
            }
        }
        else
        {
            // From master
            List<label> recv;

            IPstream::recv(recv, UPstream::masterNo());

            Pout<<"recv: " << flatOutput(recv) << endl;
        }

        UPstream::barrier(UPstream::worldComm);
    }


    // SortList
    {
        List<scalar> list1(20);

        Random rnd(1234);

        for (scalar& val : list1)
        {
            val = 100 * rnd.sample01<scalar>();
        }

        // Pick out 1/2 the values and make the negative
        for (label i=0; i < list1.size()/2; ++i)
        {
            label pos = rnd.position(0, list1.size()-1);
            list1[pos] = -list1[pos];
        }

        Info<< nl << "Random list: " << flatOutput(list1) << nl;

        SortList<scalar> sorter1(list1);

        Info<< nl << "Sorted: " << flatOutput(sorter1) << nl;
        Info<< nl << "Sort indices: " << flatOutput(sorter1.indices()) << nl;

        // Reverses indices
        sorter1.reverse();

        Info<< nl << "Reverse: " << flatOutput(sorter1.indices()) << nl;

        sorter1.reset();

        Info<< nl << "Reset indices: " << flatOutput(sorter1.indices()) << nl;

        sorter1.reverseSort();

        Info<< nl << "Reverse indices: " << flatOutput(sorter1.indices()) << nl;

        Info<< nl << "Reverse sorted: " << flatOutput(sorter1) << nl;

        sorter1.sort(std::greater<scalar>());

        SortList<scalar> sorter2(list1, std::greater<scalar>());

        sorter2.reverse();

        Info<< "sorted: ";
        for (const auto& val : sorter2)
        {
            Info<< ' ' << val;
        }
        Info<< nl;


        sorter2.sort([](scalar a, scalar b) { return mag(a) < mag(b); });

        Info<< nl << "Mag sorted: " << flatOutput(sorter2) << nl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
