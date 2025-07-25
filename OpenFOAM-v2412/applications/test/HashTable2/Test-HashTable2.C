/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Miscellaneous tests for HashTable

\*---------------------------------------------------------------------------*/

#include "HashTable.H"
#include "HashPtrTable.H"
#include "Map.H"
#include "FlatOutput.H"

using namespace Foam;


// Simple wrapper for testing purposes
class Label
{
    label data_;

public:

    Label()
    :
        data_(0)
    {}

    Label(label val)
    :
        data_(val)
    {}

    ~Label()
    {
        Info<<"delete label: " << data_ << endl;
    }

    // Some arbitrary non-const method (for testing)
    label increment()
    {
        return ++data_;
    }

    // Some arbitrary method (for testing)
    std::string info() const
    {
        return "boxed label=" + std::to_string(data_);
    }

    friend Ostream& operator<<(Ostream& os, const Label& val)
    {
        os  << val.data_;
        return os;
    }
};


template<class T, class Key, class Hash>
void printTable(const HashPtrTable<T, Key, Hash>& table)
{
    Info<< table.size() << '(' << nl;

    for (const auto& key : table.sortedToc())
    {
        const auto iter = table.find(key);

        Info
            << "    " << iter.key() << " (" << Foam::name(&(iter.key()))
            << ") : ";

        if (iter.val())
        {
            Info<< *(iter.val());
        }
        else
        {
            Info<< "nullptr";
        }

        Info<< " (" << Foam::name(iter.val()) << ")" << nl;;
    }

    Info<< ')' << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    HashTable<Label, Foam::string> table0
    ({
        {"abc", 123},
        {"kjhk", 10},
        {"kjhk2", 12}
    });

    Info<< "table0: " << table0 << nl
        << "toc: " << flatOutput(table0.toc()) << nl;

    HashTable<label, label, Hash<label>> table2
    (
        // From key/value pairs
        labelList({3, 5, 7}),
        labelList({10, 12, 16})
    );

    Info<< "table2: " << table2 << nl
        << "toc: " << flatOutput(table2.toc()) << nl;

    Map<label> table3(1);
    table3.transfer(table2);

    Info<< "table2: " << table2 << nl
        << "toc: " << flatOutput(table2.toc()) << nl;

    Info<< "table3: " << table3 << nl
        << "toc: " << flatOutput(table3.toc()) << nl;

    Map<label> table4(std::move(table3));

    Info<< "table3: " << table3 << nl
        << "toc: " << flatOutput(table3.toc()) << nl;

    Info<< "table4: " << table4 << nl
        << "toc: " << flatOutput(table4.toc()) << nl;

    {
        HashTable<Label, Foam::string> table1(0);
        table1.insert("abc", 5);
        table1.insert("def", 10);
        table1.insert("ghi", 15);
        table1.insert("jkl", 20);

        Info<< nl << "Table toc: " << flatOutput(table1.sortedToc()) << nl;

        for (const word k : { "abc" })
        {
            const auto iter = table1.cfind(k);

            if (iter.good())
            {
                Info<< "have " << k << nl
                    << "    info: " << (*iter).info() << nl
                    // Good: does not compile
                    // << "    info: " << iter->info() << nl
                    ;
            }

            auto iter2 = table1.find(k);

            if (iter2.good())
            {
                Info<< "have " << k << nl
                    << "    incr: " << (*iter2).increment() << nl
                    // Good: does not compile
                    // << "    incr: " << iter2->increment() << nl
                    ;
            }
        }

        Info<< nl
            << "Table0: " << flatOutput(table0.sortedToc()) << nl
            << "Table1: " << flatOutput(table1.sortedToc()) << nl;

        table1.merge(table0);

        Info<< "merged 0: " << flatOutput(table0.sortedToc()) << nl
            << "merged 1: " << flatOutput(table1.sortedToc()) << nl;
    }

    {
        HashPtrTable<Label> ptable0(0);
        ptable0.emplace("abc", 123),
        ptable0.emplace("kjhk", 10);
        ptable0.emplace("kjhk2", 12);

        HashPtrTable<Label> ptable1(0);
        ptable1.insert("abc", autoPtr<Label>::New(5));
        ptable1.emplace("def", 10);
        ptable1.emplace("ghi", 15);
        ptable1.emplace("jkl", 20);

        Info<< nl << "PtrTable toc: " << flatOutput(ptable1.sortedToc()) << nl;

        for (const word k : { "abc" })
        {
            const auto iter = ptable1.cfind(k);

            // Note: increment() changes contents of pointers,
            // not the pointers themselves.
            if (iter.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << name(*iter) << nl
                    << "    info: " << (*iter)->info() << nl
                    << "    incr: " << (*iter)->increment() << nl
                    ;
            }

            auto iter2 = ptable1.find(k);

            if (iter2.good())
            {
                Info<< "have " << k << nl
                    << "    incr: " << (*iter2)->increment() << nl
                    ;
            }
        }


        // Attempt the same again

        HashTable<Label*> tableView1;
        HashTable<const Label*> tableView2;

        forAllConstIters(ptable1, iter)
        {
            tableView1.insert(iter.key(), iter.val());
            tableView2.insert(iter.key(), iter.val());
        }

        Info<< nl << "Table<pointer> toc: "
            << flatOutput(tableView1.toc()) << nl;

        for (const word k : { "abc" })
        {
            const auto iter1 = tableView1.cfind(k);

            // Note that increment changes contents of the pointers
            // not the table
            if (iter1.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << name(*iter1) << nl
                    << "    info: " << (*iter1)->info() << nl
                    << "    incr: " << (*iter1)->increment() << nl
                    ;
            }

            auto iter2 = tableView2.cfind(k);
            if (iter2.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << name(*iter2) << nl
                    << "    info: " << (*iter2)->info() << nl
                    // Good: does not compile
                    // << "    incr: " << iter2->increment() << nl
                    ;
            }

            auto iter3 = tableView2.find(k);
            if (iter3.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << name(*iter3) << nl
                    << "    info: " << (*iter3)->info() << nl
                    // Good: does not compile
                    // << "    incr: " << iter3->increment() << nl
                    ;
            }
        }

        Info<< nl
            << "Table0" << nl;
        printTable(ptable0);

        Info<< "Table1" << nl;
        printTable(ptable1);

        ptable1.merge(ptable0);

        Info<< "merged 0" << nl;
        printTable(ptable0);
        Info<< "merged 1" << nl;
        printTable(ptable1);

        Info<< nl << "Ending scope" << nl;
    }

    {
        Info<< nl << "Table<labelList> copy/move/emplace insertion" << nl;

        HashTable<labelList> ltable1(0);
        ltable1.insert("abc", identity(2));
        ltable1.insert("def", identity(3));
        ltable1.insert("ghi", identity(4));
        ltable1.emplace("jkl", label(10), -35);
        ltable1.emplace("mno");
        ltable1.emplace("def", label(2), -2);      // no overwrite
        ltable1.emplace_set("ghi", label(2), -2);  // overwrite

        labelList list1(identity(4, -4));

        Info<< "move insert " << list1 << nl;

        ltable1.insert("pqr", std::move(list1));

        Info<< "after insert " << list1 << nl;

        Info<< nl << "HashTable<labelList>: "
            << ltable1 << nl;


        // Use '->' dereferencing
        const auto iter = ltable1.cfind("ghi");

        if (iter)
        {
            Info<< "got with " << (*iter).size() << nl;
        }

        Info<< nl
            << "Test (DIY) insert_or_assign" << nl;

        label nKeys = 0;

        for (const auto& key : { "abc", "foo", "mno", "xyz" })
        {
            Info<< key;
            if (ltable1.contains(key))
            {
                Info<< " : " << ltable1[key];
            }
            else
            {
                Info<< " : n/a";
            }

            /// ltable1.insert_or_assign(key, identity(++nKeys));
            ltable1(key) = identity(++nKeys);

            Info<< "  --> " << ltable1[key] << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
