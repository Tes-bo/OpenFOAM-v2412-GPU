/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "HashTable.H"
#include "List.H"
#include "DynamicList.H"
#include "FlatOutput.H"
#include "IOstreams.H"
#include "SpanStream.H"
#include "ListOps.H"
#include "stringListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    HashTable<scalar> table1
    {
        {"aaa", 1.0},
        {"aba", 2.0},
        {"aca", 3.0},
        {"ada", 4.0},
        {"aeq", 5.0},
        {"aaw", 6.0},
        {"abs", 7.0},
        {"acr", 8.0},
        {"adx", 9.0},
        {"aec", 10.0}
    };

    // Info<< "\ntable1: " << table1<< endl;

    // Erase by key
    table1.erase("aaw");

    // Erase by iterator
    HashTable<scalar>::iterator iter = table1.find("abs");
    table1.erase(iter);

    Info<< "\ntable1 toc: " << table1.toc() << endl;
    Info<< "\ntable1 sortedToc: " << table1.sortedToc() << endl;
    table1.printInfo(Info)
        << "table1 [" << table1.size() << "] " << endl;
    forAllConstIters(table1, iter)
    {
        Info<< iter.key() << " => " << iter.val() << nl;
    }


    Info<< "\ntable1 csorted() :" << endl;
    for (const auto& iter : table1.csorted())
    {
        Info<< "    " << iter.key() << " => " << iter.val() << nl;
    }
    Info<< endl;


    table1.set("acr", 108);
    table1.set("adx", 109);
    table1.set("aec", 100);
    table1("aaw") -= 1000;
    table1("aeq") += 1000;

    Info<< "\noverwrote some values table1: " << table1 << endl;


    // Test writable sorted access
    for (auto& iter : table1.sorted())
    {
        // Should not compile:  iter.key() = "illegal";
        iter.val() *= 2;
    }

    Info<< "\nInplace modified - via sorted() access :" << endl;
    for (const auto& iter : table1.csorted())
    {
        Info<< "    " << iter.key() << " => " << iter.val() << nl;
    }
    Info<< endl;


    Info<< "\ntest find:" << endl;
    Info<< table1.find("aaa")() << nl
        << table1.find("aba")() << nl
        << table1.find("aca")() << nl
        << table1.find("ada")() << nl
        << table1.find("aeq")() << nl
        << table1.find("acr")() << nl
        << table1.find("adx")() << nl
        << table1.find("aec")() << nl
        << table1["aaa"] << nl;

    {
        OCharStream os;
        os  << table1;
        HashTable<scalar> readTable(ISpanStream(os.view())());

        Info<< "Istream constructor:" << readTable << endl;
    }


    HashTable<scalar> table2(table1); // Copy
    HashTable<scalar> table3(std::move(table1)); // Move

    Info<< nl
        << "copy table1 -> table2" << nl
        << "move table1 -> table3" << nl;

    Info<< "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl
        << "\ntable3" << table3 << nl;

    Info<< "\nerase table2 by iterator" << nl;
    forAllIters(table2, iter)
    {
        Info<< "erasing " << iter.key() << " => " << iter.val() << " ... ";
        table2.erase(iter);
        Info<< "erased" << endl;
    }

    Info<< "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl
        << "\ntable3" << table3 << nl;

    table3.resize(1);
    Info<< "\nresize(1) table3" << nl;
    table3.printInfo(Info)
        << table3 << nl;

    table3.resize(10000);
    Info<< "\nresize(10000) table3" << nl;
    table3.printInfo(Info)
        << table3 << nl;

    HashTable<scalar> table4;

    table4 = table3;
    Info<< "\ncopy table3 -> table4 " << table4 << nl;

    Info<< "\nclear table4 ... ";
    table4.clear();
    Info<< "[" << table4.size() << "] " << table4 << nl;

    table1 = table3;
    Info<< "\ncopy table3 -> table1 (previously transferred)" << table1 << nl;

    Info<< "test table1 == table3 : " << (table1 == table3) << nl;
    table1.erase(table1.begin());
    Info<< "removed an element - test table1 != table3 : "
        << (table1 != table3) << nl;

    // Insert a few things into table2
    table2.set("ada", 14.0);
    table2.set("aeq", 15.0);
    table2.set("aaw", 16.0);
    table2.set("abs", 17.0);
    table2.set("adx", 20.0);

    Info<< "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl;

    label nErased = table1.erase(table2);

    Info<< "\nerase table2 keys from table1 (removed "
        << nErased << " elements)" << nl
        << "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl;


    Info<< "\ntable3" << table2
        << "\nclearStorage table2 ... ";
    table2.clearStorage();
    Info<< table2 << nl;

    table1 =
    {
        {"abc", 3.0},
        {"def", 6.0},
        {"acr", 8.0},
        {"aec", 10.0}
    };

    Info<< "\ntable1" << table1 << nl;

    Info<< "\nrange-for(table1) - returns values" << nl;
    for (const auto& it : table1)
    {
        Info<< "val:" << it << nl;
    }

    Info<< "\nrange-for(table1.keys()) - returns keys" << nl;
    for (const auto& k : table1.keys())
    {
        Info<< "key:" << k << nl;
    }

    // These do not yet work. Issues resolving the distance.
    //
    //  List<scalar> table1vals(table1.begin(), table1.end());

    {
        Info<<"distance/size: "
            << std::distance(table1.begin(), table1.end())
            << "/" << table1.size()
            << " and "
            << std::distance(table1.keys().begin(), table1.keys().end())
            << "/" << table1.keys().size()
            << nl;

        List<word> sortKeys
        (
            ListOps::create<word>
            (
                table1.keys().begin(),
                table1.keys().end(),
                identityOp{}
            )
        );
        sort(sortKeys);
        Info<<"sortKeys: " << flatOutput(sortKeys) << nl;
    }

    Info<< "\nFrom table1: " << flatOutput(table1.sortedToc()) << nl
        << "retain keys: " << flatOutput(table3.sortedToc()) << nl;

    table1.retain(table3);
    Info<< "-> " << flatOutput(table1.sortedToc()) << nl;

    Info<< "Lookup non-existent" << nl;

    Info<< table1.lookup("missing-const", 1.2345e+6)
        << "  // const-access" << nl;

    Info<< table1("missing-inadvertent", 3.14159)
        << "  // (inadvertent?) non-const access"  << nl;

    Info<< table1("missing-autovivify")
        << "  // Known auto-vivification (non-const access)" << nl;

    Info<<"\ntable1: " << table1 << endl;

    // Start again
    HashTable<scalar> table1start
    {
        {"aaa", 1.0},
        {"aba", 2.0},
        {"a_ca", 3.0},
        {"ada", 4.0},
        {"aeq_", 5.0},
        {"aaw", 6.0},
        {"abs", 7.0},
        {"a_cr", 8.0},
        {"adx", 9.0},
        {"ae_c", 10.0}
    };

    table1 = table1start;
    Info<< "\ntable has keys: "
        << flatOutput(table1.sortedToc()) << nl;

    wordRe matcher(".*_.*", wordRe::REGEX);
    table1.filterKeys
    (
        [&matcher](const word& k){ return matcher.match(k); }
    );
    Info<< "retain things matching " << matcher << " => "
        << flatOutput(table1.sortedToc()) << nl;

    table1 = table1start;
    table1.filterKeys
    (
        [&matcher](const word& k){ return matcher.match(k); },
        true
    );

    Info<< "prune things matching " << matcher << " => "
        << flatOutput(table1.sortedToc()) << nl;

    // Same, without a lambda
    table1 = table1start;
    table1.filterKeys(matcher, true);

    Info<< "prune things matching " << matcher << " => "
        << flatOutput(table1.sortedToc()) << nl;


    // Same idea, but inverted logic inside the lambda
    table1 = table1start;
    table1.filterKeys
    (
        [&matcher](const word& k){ return !matcher.match(k); },
        true
    );

    Info<< "prune things matching " << matcher << " => "
        << flatOutput(table1.sortedToc()) << nl;


    table1 = table1start;
    Info<< "\ntable:" << table1 << nl;

    table1.filterValues
    (
        [](const scalar& v){ return (v >= 5); }
    );

    Info<< "\ntable with values >= 5:" << table1 << nl;

    table1 = table1start;
    Info<< "\ntable:" << table1 << nl;

    table1.filterEntries
    (
        [&matcher](const word& k, const scalar& v)
        {
            return matcher(k) && (v >= 5);
        }
    );

    Info<< "\ntable with values >= 5 and matching " << matcher
        << table1 << nl;


    table1 = table1start;
    Info<< "\ntable:" << table1 << nl;
    Info<< "has "
        << table1.countValues([](const scalar& v) { return v >= 7; })
        << " values >= 7 with these keys: "
        << table1.tocValues([](const scalar& v) { return v >= 7; })
        << nl;


    // Start again with new value
    table2.set("ada", 14.0);
    table2.set("aeq", 15.0);
    table2.set("aaw", 16.0);

    Info<< nl << "input values" << nl;
    Info<<"table1 =  " << table1 << nl <<"table2 =  " << table2 << nl;

    Info<<"std::swap function" << nl;
    std::swap(table1, table2);
    Info<<"table1 =  " << table1 << nl <<"table2 =  " << table2 << nl;

    Info<<"swap method" << nl;
    table1.swap(table2);
    Info<<"table1 =  " << table1 << nl <<"table2 =  " << table2 << nl;

    Info<<"transfer" << nl;
    table1.transfer(table2);
    Info<<"table1 =  " << table1 << nl <<"table2 =  " << table2 << nl;

    Info<<"move assign" << nl;
    table2 = std::move(table1);
    Info<<"table1 =  " << table1 << nl <<"table2 =  " << table2 << nl;

    Info<<"move construct" << nl;
    HashTable<scalar> table1b(std::move(table2));
    Info<<"table1 =  " << table1b << nl <<"table2 =  " << table2 << nl;

    table1b.set("more", 14.0);
    table1b.set("less", 15.0);
    table1b.set("other", 16.0);

    Info<<"Test a += b " << nl;
    Info<<"a = " << flatOutput(table1b.sortedToc()) << nl
        <<"b = " << flatOutput(table3.sortedToc()) << nl;

    Info<<"=> " << (table1b += table3) << nl;

    Info<< "\nDone\n";

    return 0;
}


// ************************************************************************* //
