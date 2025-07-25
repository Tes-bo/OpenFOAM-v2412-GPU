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

Application
    Test-List

Description
    Simple tests and examples of use of List

See also
    Foam::List

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "wordRes.H"

#include "IOstreams.H"
#include "SpanStream.H"
#include "scalar.H"
#include "vector.H"

#include "labelRange.H"
#include "scalarList.H"
#include "HashOps.H"
#include "ListOps.H"
#include "IndirectList.H"
#include "SubList.H"
#include "SliceList.H"
#include "ListPolicy.H"

#include <list>
#include <numeric>
#include <functional>

// see issue #2083
#undef Foam_constructList_from_iterators

namespace Foam
{

// Verify inheritance
class MyStrings
:
    public List<string>
{
public:

    using List<string>::List;
};



namespace Detail
{
namespace ListPolicy
{

// Override on a per-type basis
template<> struct short_length<short> : std::integral_constant<int,20> {};

} // End namespace ListPolicy
} // End namespace Detail
} // End namespace Foam



using namespace Foam;

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


void printMyString(const UList<string>& lst)
{
    MyStrings slist2(lst);

    Info<<slist2 << nl;
}


template<class T>
Ostream& printListOutputType(const char* what)
{
    Info<< what
        << " (contiguous="
        << is_contiguous<T>::value << " no_linebreak="
        << Detail::ListPolicy::no_linebreak<T>::value
        << " short_length="
        << Detail::ListPolicy::short_length<T>::value << ')';

    return Info;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();

    argList::addOption("reList", "reList");
    argList::addOption("wordList", "wordList");
    argList::addOption("stringList", "stringList");
    argList::addOption("float", "xx");
    argList::addBoolOption("create", "Test ListOps::create functionality");
    argList::addBoolOption("ListList", "Test list of list functionality");
    argList::addBoolOption("flag");

    #include "setRootCase.H"

    {
        List<label> ident(15);
        Foam::identity(ident, 0);

        Info<<"Ident:";
        forAllConstIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"reverse:";
        forAllReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"const reverse:";
        forAllConstReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<< "data:" << Foam::name(ident.cdata())
            << " size:" << ident.size() << nl;


        Info<< "resize_unsafe(10)" << nl;
        ident.resize_unsafe(10);

        Info<< "data:" << Foam::name(ident.cdata())
            << " size:" << ident.size() << nl;
    }

    if (false)
    {
        labelList intlist(ICharStream("(0 1 2)")());
        Info<<"construct from Istream: " << intlist << endl;

        ICharStream("(3 4 5)")() >> static_cast<labelUList&>(intlist);
        Info<<"is >>: " << intlist << endl;

        ICharStream("(6 7 8)")() >> intlist;
        Info<<"is >>: " << intlist << endl;
    }

    List<vector> list1(ICharStream("1 ((0 1 2))")());
    Info<< "list1: " << list1 << endl;

    List<vector> list2
    {
        vector(0, 1, 2),
        vector(3, 4, 5),
        vector(6, 7, 8),
        vector(0, 1, 2),
        vector(3, 4, 5),
        vector(6, 7, 8),
    };
    Info<< "list2: " << list2 << endl;

    Info<< "forAllConstIters(list2): ";
    forAllConstIters(list2, iter) { Info<< " " << *iter; }
    Info<< endl;

    Info<< "forAllConstReverseIters(list2): ";
    forAllConstReverseIters(list2, iter) { Info<< " " << *iter; }
    Info<< endl;

    Info<< "forAllConstIters(list2): ";
    forAllIters(list2, iter) { *iter *= 2; Info<< " " << *iter; }
    Info<< endl;

    Info<< "forAllReverseIters(list2): ";
    forAllReverseIters(list2, iter) { *iter *= 0.5; Info<< " " << *iter; }
    Info<< endl;

    list1.push_back(list2);
    Info<< "list1.push_back(list2): " << list1 << endl;

    for (const vector& val : { vector(3, 4, 5), vector(10,11, 12)} )
    {
        testFind(val, list2);
    }

    list2.setSize(10, vector(1, 2, 3));
    Info<< "list2: " << list2 << endl;

    List<vector> list3(std::move(list2));
    Info<< "Move construct" << endl;
    Info<< "list2: " << list2 << nl
        << "list3: " << list3 << endl;

    List<vector> list4
    {
        vector(0, 1, 2),
        vector(3, 4, 5),
        vector(6, 7, 8)
    };
    Info<< "list4: " << list4 << endl;

    List<vector> list5
    {
        {5, 3, 1},
        {10, 2, 2},
        {8, 1, 0}
    };
    Info<< "list5: " << list5 << endl;
    list5 =
    {
        {8, 1, 0},
        {5, 3, 1},
        {10, 2, 2}

    };
    Info<< "list5: " << list5 << endl;

    list4.swap(list5);
    Info<< "Swapped via the swap() method" << endl;
    Info<< "list4: " << list4 << nl
        << "list5: " << list5 << endl;

    #ifdef Foam_constructList_from_iterators
    List<vector> list6(list4.begin(), list4.end());
    Info<< "list6: " << list6 << endl;
    #else
    Info<< "NOTE: no construction from two iterators" << endl;
    #endif

    // Subset
    const labelList map{0, 2};
    List<vector> subList3(list3, map);
    Info<< "Elements " << map << " out of " << list3
        << " => " << subList3 << endl;

    // test flattened output
    {
        Info<< nl;

        labelList longLabelList = identity(15);

        // This will not work:
        // scalarList slist = identity(15);
        //
        // More writing, but does work:
        #ifdef Foam_constructList_from_iterators
        scalarList slist(labelRange().begin(), labelRange(15).end());

        Info<<"scalar identity:" << flatOutput(slist) << endl;
        #else
        Info<<"No iterator means of creating a scalar identity list" << endl;
        #endif

        printListOutputType<label>("labels") << nl;

        Info<< "normal: " << longLabelList << nl;
        Info<< "flatOutput: " << flatOutput(longLabelList) << nl;
        // Info<< "flatOutput(14): " << flatOutput(longLabelList, 14) << nl;

        auto shrtList = ListOps::create<short>
        (
            longLabelList,
            [](const label& val){ return val; }
        );

        printListOutputType<short>("short") << nl;
        Info<< "normal: " << shrtList << nl;


        stringList longStringList(12);
        forAll(longStringList, i)
        {
            longStringList[i].resize(3, 'a' + i);
        }

        printListOutputType<string>("string") << nl;

        Info<< "normal: " << longStringList << nl;
        Info<< "flatOutput: " << flatOutput(longStringList) << nl;

        auto wList = ListOps::create<word>
        (
            longStringList,
            [](const std::string& val){ return val; }
        );

        printListOutputType<word>("word") << nl;

        Info<< "normal: " << wList << nl;

        // Shorten
        longStringList.resize(8);
        wList.resize(8);

        Info<< "Test shorter lists" << nl;

        printListOutputType<string>("string") << nl;
        Info<< "normal: " << longStringList << nl;

        printListOutputType<word>("word") << nl;
        Info<< "normal: " << wList << nl;
    }

    // Test SubList and labelRange
    {
        Info<< nl;
        labelList longLabelList = identity(25);
        reverse(longLabelList);

        FixedList<label, 6> fixedLabelList({0,1,2,3,4,5});
        const labelList constLabelList = identity(25);

        Info<< "full-list: " << flatOutput(longLabelList) << nl;

        labelRange range1(-15, 25);
        Info<<"sub range:" << range1 << "=";
        Info<< SubList<label>(longLabelList, range1) << nl;

        {
            // A valid range
            const labelRange subset(4, 5);

            // Assign some values
            longLabelList.slice(subset) = identity(subset.size());

            Info<<"assigned identity in range:" << subset
                << "=> " << flatOutput(longLabelList) << nl;

            // Assign values in iterator range
            std::iota(longLabelList.begin(15), longLabelList.begin(50), 115);

            Info<<"assigned values in iterator range "
                << "=> " << flatOutput(longLabelList) << nl;

            labelList someList(identity(24));

            longLabelList.slice(subset) =
                SliceList<label>(someList, sliceRange(8, subset.size(), 2));

            Info<<"assigned sliced/stride in range:" << subset
                << "=> " << flatOutput(longLabelList) << nl;

            // Does not work - need a reference, not a temporary
            // Foam::reverse(longLabelList[subset]);

            {
                auto sub(longLabelList.slice(subset));
                Foam::reverse(sub);
            }

            Info<<"reversed range:" << subset
                << "=> " << flatOutput(longLabelList) << nl;
        }

        labelRange range2(7, 8);
        Info<<"sub range:" << range2 << "=";
        Info<< SubList<label>(longLabelList, range2) << nl;

        // labelRange range2(7, 8);
        Info<<"use range " << range2 << " to set value";
        SubList<label>(longLabelList, range2) = -15;
        Info<< "=> " << flatOutput(longLabelList) << nl;

        // This syntax looks even nicer:

        // GOOD: does not compile
        // > constLabelList[labelRange(23,5)] = 5;

        // Check correct overlaps
        longLabelList.slice(labelRange(-10, 12)) = 200;
        longLabelList.slice({18,3}) = 100;
        longLabelList.slice({23,3}) = 400;
        // and complete misses
        longLabelList.slice({500,50}) = 100;

        // -ve size suppressed by internal 'validateRange' = no-op
        longLabelList.slice({5,-5}) = 42;
        longLabelList.slice({21,100}) = 42;

        #if 0
        // Compiles, but is runtime abort!
        const bool oldThrowingError = FatalError.throwing(true);
        try
        {
            longLabelList.slice(labelRange(20,50)) = constLabelList;
        }
        catch (const Foam::error& err)
        {
            Info<< "Caught FatalError "
                << err << nl << endl;
        }
        FatalError.throwing(oldThrowingError);
        #endif

        //Good: does not compile
        // longLabelList[labelRange(20,50)] = fixedLabelList;

        Info<< "updated: " << constLabelList.slice(labelRange(23,5)) << nl;
        Info<< "updated: " << flatOutput(longLabelList) << nl;

        //Nope: sort(longLabelList.slice(labelRange(18,5)));
        {
            // Instead
            auto sub = longLabelList.slice(labelRange(8));
            sort(sub);
        }
        Info<< "sub-sorted: " << flatOutput(longLabelList) << nl;

        #ifdef Foam_constructList_from_iterators
        // Construct from a label-range
        labelRange range(25,15);

        labelList ident(range.begin(), range.end());
        Info<<"range-list (label)=" << ident << nl;

        List<scalar> sident(range.begin(), range.end());
        Info<<"range-list (scalar)=" << sident << nl;

//        // Sub-ranges also work
//        List<scalar> sident2(range.at(3), range.at(10));
//        Info<<"subrange-list (scalar)=" << sident2 << nl;

        // VERY BAD IDEA: List<scalar> sident3(range.at(10), range.at(3));

        // This doesn't work, and don't know what it should do anyhow
        // List<vector> vident(range.begin(), range.end());
        // Info<<"range-list (vector)=" << vident << nl;

        // Even weird things like this
        List<scalar> sident4(labelRange().begin(), labelRange(8).end());
        Info<<"range-list (scalar)=" << sident4 << nl;
        #else
        Info<< "NOTE: no construction of labelList from range pair" << nl
            << "use identity(...) instead" << endl;
        #endif
    }

    wordReList reLst;
    wordList wLst;
    stringList sLst;

    scalar xxx(-1);

    if (args.found("flag"))
    {
        Info<<"-flag:" << args["flag"] << endl;
    }

    if (args.found("create"))
    {
        Info<< nl << "Test ListOps::create functionality" << nl;

        const auto labels = identity(15);
        Info<< "labels: " << flatOutput(labels) << endl;

        {
            auto scalars = ListOps::create<scalar>
            (
                labels,
                [](const label& val){ return scalar(1.5*val); }
            );
            Info<< "scalars: " << flatOutput(scalars) << endl;
        }

        {
            auto vectors = ListOps::create<vector>
            (
                labels,
                [](const label& val){ return vector(1.2*val, -1.2*val, 0); }
            );
            Info<< "vectors: " << flatOutput(vectors) << endl;
        }

        {
            auto longs = ListOps::create<long>
            (
                labels,
                [](const label& val){ return val; }
            );
            Info<< "longs: " << flatOutput(longs) << endl;
        }
        {
            auto negs = ListOps::create<label>
            (
                labels,
                std::negate<label>()
            );
            Info<< "negs: " << flatOutput(negs) << endl;
        }

        {
            auto scalars = ListOps::create<scalar>
            (
                labelRange().cbegin(),
                labelRange(15).cend(),
                [](const label& val){ return scalar(-1.125*val); }
            );
            Info<< "scalars: " << flatOutput(scalars) << endl;
        }

        // #if WM_LABEL_SIZE == 32
        {
            List<int64_t> input(10);
            Foam::identity(input, 50);

            auto output = ListOps::create<label>
            (
                input,
                labelOp<int64_t>()
            );
            Info<< "label (from int64): " << flatOutput(output) << endl;
        }
        // #elif WM_LABEL_SIZE == 64
        {
            List<int32_t> input(10);
            Foam::identity(input, 50);

            auto output = ListOps::create<label>
            (
                input,
                labelOp<int32_t>()
            );
            Info<< "label (from int32): " << flatOutput(output) << endl;
        }
        // #endif


        labelHashSet locations{ -15, 5, 10, 15, 25, 35 };
        Info<< nl << "Test for createWithValue with locations :"
            << flatOutput(locations.sortedToc()) << nl;

        {
            auto output = ListOps::createWithValue<label>
            (
                30,
                locations.toc(),  // Any order
                100,
                -1  // default value
            );
            Info<< "with labelUList: " << flatOutput(output)
                << " selector: " << flatOutput(locations.sortedToc()) << nl;
        }

        {
            auto output = ListOps::createWithValue<label>
            (
                30,
                locations,
                100,
                -1  // default value
            );
            Info<< "with labelHashSet: " << flatOutput(output)
                << " selector: " << flatOutput(locations) << nl;
        }

        {
            bitSet select = HashSetOps::bitset(locations);
            auto output = ListOps::createWithValue<label>
            (
                30,
                select,
                100,
                -1  // default value
            );
            Info<< "with bitSet: " << flatOutput(output)
                << " selector: " << flatOutput(select.toc()) << nl;
        }

        {
            List<bool> select = HashSetOps::bools(locations);

            auto output = ListOps::createWithValue<label>
            (
                30,
                select,
                100,
                -1  // default value
            );
            Info<< "with boolList: " << flatOutput(output)
                << " selector: " << flatOutput(select) << nl;
        }

        // Repeat with a shorter selector
        locations = { -15, 5, 10 };

        {
            auto output = ListOps::createWithValue<label>
            (
                30,
                locations,
                100,
                -1  // default value
            );
            Info<< "with labelHashSet: " << flatOutput(output)
                << " selector: " << flatOutput(locations) << nl;
        }

        {
            bitSet select = HashSetOps::bitset(locations);
            auto output = ListOps::createWithValue<label>
            (
                30,
                select,
                100,
                -1  // default value
            );
            Info<< "with bitSet: " << flatOutput(output)
                << " selector: " << flatOutput(HashSetOps::used(select))
                << nl;
        }

        {
            List<bool> select = HashSetOps::bools(locations);

            auto output = ListOps::createWithValue<label>
            (
                30,
                select,
                100,
                -1  // default value
            );
            Info<< "with boolList: " << flatOutput(output)
                << " selector: " << flatOutput(HashSetOps::used(select))
                << nl;
        }
    }

    if (args.found("ListList"))
    {
        {
            labelListList listlist(5, identity(5));
            Info<<"list-list with length/val:" << listlist << nl;
        }

        {
            labelListList listlist(one{}, identity(5));
            Info<<"list-list 1/val:" << listlist << nl;
        }

        {
            labelList content = identity(5);

            labelListList listlist(one{}, content);
            Info<<"list-list 1/copy val:" << listlist
                <<" - from " << content << nl;
        }

        {
            labelList content = identity(5);

            labelListList listlist(one{}, std::move(content));
            Info<<"list-list 1/move val:" << listlist
                <<" - from " << content << nl;
        }

        {
            labelListList listlist(one{}, Zero);
            Info<<"list-list 1/move val:" << listlist
                << nl;
        }
    }

    if (args.readIfPresent<scalar>("float", xxx))
    {
        Info<<"read float " << xxx << endl;
    }

    args.readListIfPresent<wordRe>("reList", reLst);
    args.readListIfPresent<word>("wordList", wLst);

    if (args.readListIfPresent<string>("stringList", sLst))
    {
        printMyString(sLst);
    }

    Info<< nl
        << "-reList:     " << flatOutput(reLst) << nl
        << "-wordList:   " << flatOutput(wLst)  << nl
        << "-stringList: " << flatOutput(sLst)  << endl;

    // Hash values
    {
        labelList list1(identity(5));
        labelList list2(identity(5));

        Info<<"hash of " << flatOutput(list1)
            << " = " << Hash<labelList>()(list1) << " or "
            << labelList::hasher()(list1) << nl;

        Info<<"hash of " << flatOutput(list2) << " = "
            << Hash<labelList>()(list2) << " or "
            << labelList::hasher()(list2) << nl;
    }

    return 0;
}

// ************************************************************************* //
