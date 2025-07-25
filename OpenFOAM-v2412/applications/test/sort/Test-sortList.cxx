/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Tests for SortList and Sortable List

\*---------------------------------------------------------------------------*/

#include "SortList.H"
#include "SortableList.H"
#include "IndirectList.H"
#include "ListOps.H"
#include "HashSet.H"
#include "stringOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:


template<class T>
void printInfo(const SortableList<T>& list)
{
    Info<< "sorted:   " << list << nl
        << "indices:  " << list.indices() << endl;
}


int main(int argc, char *argv[])
{
    {
        const labelList orig({7, 9, 7, 12, 9, 1, 12, 2, 4, 7, 4, 0});

        Info<< "original list: " << flatOutput(orig) << endl;

        auto unique = IndirectList<label>::uniq(orig);

        Info<< "unique  list: " << flatOutput(unique) << endl;
        Info<< "unique  idx : " << flatOutput(unique.addressing()) << endl;
    }


    const labelList orig({7, 9, 1, 2, 4, 7, 4, 0});

    labelList order;

    labelList list1(orig);
    sortedOrder(list1, order);

    SortableList<label> list1r(list1.size());
    list1r = list1;

    Info<< "unsorted: " << orig << nl
        << "order:    " << sortedOrder(list1) << endl;
    sort(list1);
    Info<< "sorted:   " << list1 << nl
        << "indices:  " << order << nl
        << "order:    " << sortedOrder(list1) << endl;

    list1r.reverseSort();
    Info<< "reverse ..." << nl;
    printInfo(list1r);

    list1r.partialSort(list1r.size()/2);
    Info<< "partial sorted ..." << nl;
    printInfo(list1r);

    SortableList<label> list2(orig);
    Info<< "unsorted: " << orig << nl;
    printInfo(list2);

    Info<< "shrunk:   " << list2.shrink() << endl;
    Info<< "indices:  " << list2.indices() << endl;

    // repeat by assignment
    list2 = orig;
    Info<< "unsorted: " << list2 << endl;
    list2.sort();

    printInfo(list2);

    // find unique/duplicate values
    list2 = orig;

    Info<< "unsorted: " << list2 << endl;
    uniqueOrder(list2, order);
    Info<< "unique:  " << order << endl;
    duplicateOrder(list2, order);
    Info<< "duplicate:" << order << endl;

    // sort reverse
    Info<< "unsorted: " << list2 << endl;
    list2.reverseSort();
    Info<< "rsort:    " << list2 << endl;
    Info<< "indices:  " << list2.indices() << endl;

    // transfer assignment
    {
        list1 = orig;
        list2.transfer(list1);
        Info<< "unsorted: " << list2 << endl;
        list2.sort();

        printInfo(list2);

        list1.transfer(list2);

        Info<< "plain:    " << list1 << endl;
        printInfo(list2);
    }

    // sort/duplicate/unique with identical values
    list2.setSize(8);
    list2 = 5;

    Info<< "unsorted: " << list2 << endl;

    uniqueOrder(list2, order);
    Info<< "unique:  " << order << endl;
    duplicateOrder(list2, order);
    Info<< "duplicate:" << order << endl;
    list2.sort();

    printInfo(list2);

    // with a single value
    list2.setSize(1);

    Info<< "unsorted: " << list2 << endl;
    uniqueOrder(list2, order);
    Info<< "unique:  " << order << endl;
    duplicateOrder(list2, order);
    Info<< "duplicate:" << order << endl;
    list2.sort();

    printInfo(list2);

    {
        labelList tmp(orig);

        Info<< nl << "move construct from List: " << tmp << endl;
        SortableList<label> list3(std::move(tmp));

        Info<< "input:    " << tmp << endl;
        printInfo(list3);

        list3.reverseSort();

        Info<< nl << "move construct from SortableList: " << list3 << endl;

        SortableList<label> list4(std::move(list3));
        Info<< "input:   " << list3 << endl;
        printInfo(list3);
        printInfo(list4);

        tmp = orig;

        Info<< nl << "move assign from List: " << tmp << endl;
        list3 = std::move(tmp);

        Info<< "input:    " << tmp << endl;
        printInfo(list3);

        Info<< nl << "move assign from SortableList: " << list3 << endl;
        list4 = std::move(list3);

        Info<< "input:    " << list3 << endl;
        printInfo(list3);
        printInfo(list4);

        labelList values;
        Info<< "move to flat-list: " << list4 << endl;

        values = std::move(list4);
        Info<< "input:    " << list4 << endl;
        printInfo(list4);
        Info<< "flat = " << values << endl;
    }

    // Sort strings
    {
        HashSet<string> hashed
        {
            "2.txt",
            "05.txt",
            "15.txt",
            "other.bak04",
            "other.bak1",
            "file1.txt",
            "file10.txt",
            "file2.txt",
            "file100.txt",
            "file.txt",
            "file011.txt",
            "file15.txt",
            "file0009.txt",
            "abcd.txt",

         // Some regular processor directories
            "processor0",
            "processor1",
            "processor9",
            "processor10",
            "processor11",
            "processor20",
            "processor21",
            "processor35",
            "processors",

         // Aggregate processor directories
            "processor0-31",
            "processor32-63",
            "processor64-95",
            "processor96-127",
            "processor128-159",
            "processor160-191",
            "processor192-223",
            "processor224-255",
        };

        Info<< nl << "Test string sorting" << nl << endl;

        // Using hash toc
        if (true)
        {
            Info<< "Unsorted" << hashed.toc() << endl;
            Info<< "sortedToc" << hashed.sortedToc() << endl;
            Info<< "natural"
                << hashed.sortedToc(stringOps::natural_sort()) << endl;

            Info<< "reverse natural"
                << hashed.sortedToc(stringOps::natural_sort::reverse())
                << endl;
        }

        // Normal list
        if (true)
        {
            labelList order;

            List<string> strings(hashed.toc());
            Info<< nl << "stringList:" << strings << endl;

            sort(strings);
            Info<< "normal sort:" << strings << endl;

            shuffle(strings);
            sort(strings, stringOps::natural_sort());
            Info<< "natural sort:" << strings << endl;

            shuffle(strings);
            sort(strings, stringOps::natural_sort::reverse());
            Info<< "reverse natural:" << strings << endl;

            strings = hashed.toc();

            Info<< nl << "test sorted order" << endl;
            Info<< nl << "list:" << strings << endl;

            sortedOrder(strings, order);
            Info<< "sortedOrder:" << flatOutput(order) << endl;

            shuffle(strings);
            sort(strings, stringOps::natural_sort());
            Info<< "reverse natural:" << strings << endl;

            shuffle(strings);
            sort(strings, stringOps::natural_sort::reverse());
            Info<< "reverse natural:" << strings << endl;

            sortedOrder
            (
                strings,
                order,
                stringOps::natural_sort::less<string>(strings)
            );
            Info<< "natural sortedOrder: " << flatOutput(order) << endl;
        }

        // SortableList
        if (false)
        {
            SortableList<string> sortable;
            Info<< nl << "Testing sortable list";

            // Assign to ensure list is initially unsorted
            sortable = hashed.toc();
            Info<< nl << "input:" << sortable << endl;

            sortable.sort();
            Info<< nl << "normal:" << sortable << endl;

            // This is still a bother (looks fairly ugly)
            // so not implemented for now

            ///  // Assign to ensure list is initially unsorted
            ///  sortable = hashed.toc();
            ///  sortable.sort
            ///  (
            ///      stringOps::natural_sort::less<string>(sortable)
            ///  );
            ///  Info<< nl << "natural:" << sortable << endl;

            ///  // Assign to ensure list is initially unsorted
            ///  sortable = hashed.toc();
            ///  sortable.sort
            ///  (
            ///      stringOps::natural_sort::greater<string>(sortable)
            ///  );
            ///  Info<< nl << "natural:" << sortable << endl;
        }

    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
