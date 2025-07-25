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

Application
    Test-bitSet2

Description
    Test bitSet functionality

\*---------------------------------------------------------------------------*/

#include "uLabel.H"
#include "boolList.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "ITstream.H"
#include "SpanStream.H"
#include "bitSet.H"
#include "FlatOutput.H"

// #define TEST_SFINAE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

inline Ostream& extent(const bitSet& bitset)
{
    Info<< "first: " << bitset.find_first()
        << " last: " << bitset.find_last()
        << " first_not: " << bitset.find_first_not()
        << endl;

    return Info;
}


inline Ostream& info(const bitSet& bitset)
{
    Info<< "size=" << bitset.size() << "/" << bitset.capacity()
        << " count=" << bitset.count()
        << " !count=" << bitset.count(false)
        << " all:" << bitset.all()
        << " any:" << bitset.any()
        << " none:" << bitset.none() << nl;

    return Info;
}


inline Ostream& info(const UList<bool>& bools)
{
    Info<< "size=" << bools.size()
        << " count=" << BitOps::count(bools)
        << " !count=" << BitOps::count(bools, false)
        << " all:" << BitOps::all(bools)
        << " any:" << BitOps::any(bools)
        << " none:" << BitOps::none(bools) << nl;

    return Info;
}


inline Ostream& report
(
    const bitSet& bitset,
    bool showBits = false,
    bool debugOutput = false
)
{
    info(bitset);

    Info<< "values: " << flatOutput(bitset) << nl;
    if (showBits)
    {
        bitset.printBits(Info, debugOutput) << nl;
    }

    return Info;
}


inline Ostream& report(const UList<bool>& bools)
{
    info(bools);
    return Info;
}


template<class UIntType>
std::string toString(UIntType value, char off='.', char on='1')
{
    std::string str(std::numeric_limits<UIntType>::digits, off);

    unsigned n = 0;

    // Starting from most significant bit - makes for easy reading.
    for
    (
        unsigned test = (1u << (std::numeric_limits<UIntType>::digits-1));
        test;
        test >>= 1u
    )
    {
        str[n++] = ((value & test) ? on : off);
    }

    return str;
}


inline bool compare
(
    const bitSet& bitset,
    const std::string& expected
)
{
    const List<unsigned int>& store = bitset.storage();

    std::string has;

    for (label blocki=0; blocki < bitset.nBlocks(); ++blocki)
    {
        has += toString(store[blocki]);
    }

    if (has == expected)
    {
        Info<< "pass: " << has << nl;
        return true;
    }

    Info<< "fail:   " << has << nl;
    Info<< "expect: " << expected << nl;

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<< "bitSet::null()" << nl
        << "sizeof : " << sizeof(bitSet::null()) << " bytes" << nl;

    info(bitSet::null());

    {
        bitSet emptySet;
        info(emptySet);
        extent(emptySet);

        emptySet.resize(10);
        info(emptySet);
        extent(emptySet);
    }

    Info<< nl << "Tests" << nl;

    bitSet list1(22);
    // Set every third one on
    forAll(list1, i)
    {
        list1[i] = !(i % 3);
    }
    Info<< "\nalternating bit pattern\n";
    compare(list1, "..........1..1..1..1..1..1..1..1");

    // As boolList
    {
        boolList bools = list1.values();

        Info<<"===============" << nl;
        Info<<"bools: " << flatOutput(bools) << nl;

        for (int i : { -10, 0, 8, 15, 32})
        {
            Info<< i << " is " << (bools.test(i) ? "set" : "unset")
                << " = " << bools.get(i) << nl;
        }

        Info<<"bools: " << flatOutput(bools) << nl;

        for (int i : { -10, 5, 24})
        {
            Info<< "set(" << i << ") = "
                << (bools.set(i) ? "true" : "false")
                << nl;
        }
        Info<<"bools: " << flatOutput(bools) << nl;

        for (int i : { -10, 12, 32, 150})
        {
            Info<< "unset(" << i << ") = "
                << (bools.unset(i) ? "true" : "false")
                << nl;
        }

        Info<<"bools: " << flatOutput(bools) << nl;

        #if 0
        boolList bools2(6, false);
        Info<<"other bools: " << flatOutput(bools2) << nl;

        bools2.set(4);
        Info<<"other bools: " << flatOutput(bools2) << nl;

        bools2.clear();
        bools2.set(3);
        bools2.resize(8);
        Info<<"other bools: " << flatOutput(bools2) << nl;
        #endif
        Info<<"===============" << nl;
    }

    #ifdef TEST_SFINAE
    {
        labelList labels = list1.toc();
        if (labels.test(0))
        {
            Info<<"no" << endl;
        }

        List<double*> ptrs(10, nullptr);
        if (ptrs.get(0))
        {
            Info<<"no" << endl;
        }
    }
    #endif


    list1.unset(labelRange(13, 20));  // In range

    Info<< "\nafter clear [13,..]\n";
    compare(list1, "...................1..1..1..1..1");

    report(list1, true);

    list1.unset(labelRange(40, 20));  // out of range
    Info<< "\nafter clear [40,..]\n";
    compare(list1, "...................1..1..1..1..1");

    report(list1, true);
    extent(list1);

    Info<< "iterate through:";
    for (const label idx : list1)
    {
        Info<<" " << idx;
    }
    Info<< nl;

    Info<< "\nalternating bit pattern\n";
    report(list1, true);

    bitSet list2 = ~list1;

    Info<< "\nflipped bit pattern\n";
    report(list2, true);
    extent(list2);

    Info<< "\nsparse set\n";
    {
        bitSet sparse(1000);
        sparse.set(300);

        info(sparse);
        extent(sparse);

        sparse.set(0);

        info(sparse);
        extent(sparse);
    }

    // set every other on
    forAll(list2, i)
    {
        list2[i] = !(i % 2);
    }

    Info<< "\nstarting pattern\n";
    report(list2, true);

    list2.resize(28, false);
    list2.resize(34, true);
    list2.resize(40, false);
    for (label i=0; i < 4; ++i)
    {
        list2[i] = true;
    }

    Info<< "\nresized with false, [28,34) true + 6 false, bottom 4 bits true\n";

    compare
    (
        list2,
        "1111.......1.1.1.1.1.1.1.1.11111"
        "..............................11"
    );


    report(list2, true);

    labelList list2Labels = list2.toc();

    Info<< "\noperator|\n";

    list1.printBits(Info);
    list2.printBits(Info);
    Info<< "==\n";
    (list1 | list2).printBits(Info);

    Info<< "\noperator& : does trim\n";
    report((list1 & list2), true);

    Info<< "\noperator^\n";
    report((list1 ^ list2), true);

    Info<< "\noperator|=\n";
    {
        bitSet a; a.set(labelRange(15, 8));
        bitSet b; b.set(labelRange(30, 8));

        report(a, false);
        report(b, false);
        report((a |= b), true);
    }
    {
        Info<< " other test\n";

        bitSet list3 = list1;
        report((list3 |= list2), true);
    }

    Info<< "\noperator&=\n";
    {
        bitSet a; a.set(labelRange(30, 8));
        bitSet b; b.set(labelRange(0, 8));

        report(a, false);
        report(b, false);
        Info<< "  (a & b)\n";
        report((a & b), true);
        Info<< "  (b & a)\n";
        report((b & a), true);
    }
    {
        Info<< " other test\n";

        bitSet list3 = list1;
        report((list3 &= list2), true);
    }

    Info<< "\noperator^=\n";
    {
        bitSet list3 = list1;
        report((list3 ^= list2), true);
    }

    Info<< "\noperator-=\n";
    {
        bitSet list3 = list1;
        report((list3 -= list2), true);
    }

    bitSet list4
    (
        ITstream
        (
            "(1 n 1 n 1 n 1 1 off 0 0 f f 0 y yes y true y false on t)"
        )()
    );

    Info<< "\ntest Istream constructor\n";

    report(list4, true);

    Info<< "\nclear/assign from labelList\n";
    list4.clear();
    list4.set(labelList{0, 1, 2, 3, 12, 13, 14, 19, 20, 21});

    report(list4, true);

    // Not yet:
    // bitSet list5{0, 1, 2, 3, 12, 13, 14, 19, 20, 21};
    // list5.printInfo(Info, true);
    // Info<< list5 << " indices: " << list5.toc() << nl;

    Info<< "\nassign from indices\n";
    list4.readList
    (
        ICharStream
        (
            "{0 1 2 3 12 13 14 19 20 21}"
        )()
    );

    report(list4, true);
    compare(list4, "..........111....111........1111");

    list4.set(labelRange(28, 6));  // extends size

    Info<<"extended\n";
    compare
    (
        list4,
        "1111......111....111........1111"
        "..............................11"
    );

    list4.set(labelRange(40, 6));  // extends size
    Info<<"extended\n";
    compare
    (
        list4,
        "1111......111....111........1111"
        "..................111111......11"
    );

    list4.unset(labelRange(14, 19));
    Info<<"cleared [14,33)\n";
    compare
    (
        list4,
        "..................11........1111"
        "..................111111......1."
    );


    // Test begin vs find_first etc
    {
        // Clear some values
        list4.unset(labelRange(0, 10));

        Info<< nl
            << "Test first vs begin" << nl
            << "    values:" << flatOutput(list4.toc()) << nl
            << "    first:" << list4.find_first() << nl
            << "    begin:" << *(list4.begin()) << nl
            << "    begin(0):" << *(list4.begin(0)) << nl;

        // With offset

        Info<< nl
            << "Test first vs begin" << nl
            << "    next(35):" << list4.find_next(35)
            << "    begin(35):" << *(list4.begin(35)) << nl
            << "    next(40):" << list4.find_next(40)
            << "    begin(40):" << *(list4.begin(40)) << nl
            << nl;
    }


    // Construct from labelUList, labelUIndList
    {
        DynamicList<label> indices({10, 50, 300});

        Info<< "set: " << flatOutput(indices) << endl;

        bitSet bools1(indices);

        Info<< "used: " << bools1.size() << " "
            << flatOutput(bools1.toc()) << endl;
    }


    // Check bitSet vs boolList
    {
        boolList bools(list4.values());

        Info<< nl << "Check BitOps on boolList" << nl << nl;

        Info<<"bitSet ";
        report(list4);

        list4.shrink();
        Info<<"shrunk ";
        report(list4);

        Info<< nl;

        Info<<"bools ";
        report(bools);

        bools = false;
        Info<<"bools (all unset) ";
        report(bools);

        bools = true;
        Info<<"bools (all set) ";
        report(bools);
    }

    Info<< "\nDone" << nl << endl;
    return 0;
}


// ************************************************************************* //
