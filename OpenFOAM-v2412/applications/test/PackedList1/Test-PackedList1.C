/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "uLabel.H"
#include "labelRange.H"
#include "bitSet.H"
#include "FlatOutput.H"
#include "IOstreams.H"

using namespace Foam;

template<unsigned Width>
inline Ostream& report
(
    const PackedList<Width>& list,
    bool showBits = false,
    bool debugOutput = false
)
{
    Info<< list.info();
    if (showBits)
    {
        list.printBits(Info, debugOutput) << nl;
    }

    return Info;
}


inline Ostream& report
(
    const bitSet& bitset,
    bool showBits = false,
    bool debugOutput = false
)
{
    Info<< bitset.info();
    Info<< "all:" << bitset.all()
        << " any:" << bitset.any()
        << " none:" << bitset.none() << nl;

    if (showBits)
    {
        bitset.printBits(Info, debugOutput) << nl;
    }

    return Info;
}

// BitOps::printBits((Info<< list1.info()), true) << nl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<< "\ntest allocation with value\n";
    PackedList<3> list1(5,1);
    report(list1);

    Info<< "\ntest assign uniform value\n";
    list1 = 3;
    report(list1);

    Info<< "\ntest assign uniform value (with overflow)\n";
    list1 = -1;
    report(list1);

    Info<< "\ntest zero\n";
    list1 = 0;
    report(list1);

    Info<< "\ntest set() with default argument (max_value)\n";
    list1.set(1);
    list1.set(3);
    report(list1);

    Info<< "\ntest unset() with in-range and out-of-range\n";
    list1.unset(3);
    list1.unset(100000);
    report(list1);

    Info<< "\ntest assign between references\n";
    list1[2] = 3;
    list1[4] = list1[2];
    report(list1);

    Info<< "\nset auto-vivify entries\n";
    list1[1] = 2;
    list1.set(8, 2);
    list1.set(10, 2);
    list1.set(14, 2);
    report(list1);

    Info<< "values()    : " << flatOutput(list1.unpack<char>()) << nl
        << "values(5,8) : " << flatOutput(list1.unpack<char>(labelRange(5,8)))
        << nl;

    {
        labelList locations({-5, -2, 2, 1, 8});

        Info<< "values at " << flatOutput(locations) << " = "
            << flatOutput(list1.unpack<char>(locations))
            << nl;
    }

    Info<< "\ntest operator== between references\n";
    if (list1[1] == list1[8])
    {
        Info<< "[1] == [8] (expected)\n";
    }
    else
    {
        Info<< "[1] != [8] (unexpected)\n";
    }

    if (list1[0] != list1[1])
    {
        Info<< "[0] != [1] (expected)\n";
    }
    else
    {
        Info<< "[0] == [1] (unexpected)\n";
    }

    {
        const PackedList<3>& constLst = list1;
        Info<< "\ntest operator[] const with out-of-range index\n";

        report(constLst);

        if (constLst[20])
        {
            Info<< "[20] is true (unexpected)\n";
        }
        else
        {
            Info<< "[20] is false (expected) list size should be unchanged "
                << "(const)\n";
        }
        report(constLst);

        Info<< "\ntest operator[] non-const with out-of-range index\n";

        // Expect failure
        const bool oldThrowingError = FatalError.throwing(true);

        try
        {
            if (list1[20])
            {
                Info<< "[20] is true (unexpected)\n";
            }
            else
            {
                Info<< "[20] is false (expected) but list was resized?? "
                    << "(non-const)\n";
            }
        }
        catch (const Foam::error& err)
        {
            Info<< "Failed (expected) " << err << nl << endl;
        }

        FatalError.throwing(oldThrowingError);
        report(list1);
    }

    {
        Info<< "\ntest operator[] with out-of-range index\n";

        // Expect failure
        const bool oldThrowingError = FatalError.throwing(true);

        try
        {
            if (!list1[20])
            {
                Info<< "[20] is false, as expected for const-access\n";
            }
        }
        catch (const Foam::error& err)
        {
            Info<< "Failed (expected) " << err << nl << endl;
        }

        FatalError.throwing(oldThrowingError);
        report(list1);
    }

    Info<< "\ntest resize with value (without reallocation)\n";
    list1.resize(8, list1.max_value);
    report(list1);

    Info<< "\ntest set() function\n";
    list1.set(1, 5);
    report(list1);

    Info<< "\ntest assign bool\n";
    list1 = false;
    report(list1);

    Info<< "\ntest assign bool\n";
    list1 = true;
    report(list1);

    Info<< "\ntest resize without value (with reallocation)\n";
    list1.resize(12);
    report(list1);

    Info<< "\ntest resize with value (with reallocation)\n";
    list1.resize(25, list1.max_value);
    report(list1);

    Info<< "\ntest resize smaller (should not touch allocation)\n";
    list1.resize(8);
    report(list1);

    Info<< "\ntest push_back() operation\n";
    list1.push_back(2);
    list1.push_back(3);
    list1.push_back(4);
    report(list1);

    Info<< "\ntest reserve() operation\n";
    list1.reserve(32);
    report(list1);

    Info<< "\ntest shrink() operation\n";
    list1.shrink();
    report(list1);

    Info<< "\ntest setCapacity() operation\n";
    list1.setCapacity(15);
    report(list1);

    Info<< "\ntest setCapacity() operation\n";
    list1.setCapacity(100);
    report(list1);

    // Expect failure
    {
        const bool oldThrowingError = FatalError.throwing(true);

        Info<< "\ntest operator[] assignment with auto-vivify\n";

        try
        {
            list1[16] = 5;
            list1[36] = list1.max_value;
        }
        catch (const Foam::error& err)
        {
            Info<< "Failed (expected) " << err << nl << endl;

            Info<< "Using set(...) instead" << nl;

            list1.set(36, list1.max_value);
        }

        FatalError.throwing(oldThrowingError);
        report(list1);
    }


    Info<< "\ntest setCapacity smaller\n";
    list1.setCapacity(24);
    report(list1);

    Info<< "\ntest resize much larger\n";
    list1.resize(150);
    report(list1);

    Info<< "\ntest trim\n";
    list1.trim();
    report(list1);

    // Add in some misc values
    list1.set(31, 1);
    list1.set(32, 2);
    list1.set(33, 3);

    Info<< "\ntest get() method\n";
    Info<< "get(10):" << list1.get(10) << " and list[10]:" << list1[10] << "\n";
    report(list1);

    Info<< "\ntest set() auto-vivify\n";
    Info<< "size:" << list1.size() << "\n";

    Info<< "list[45]:" << list1.get(45) << "\n";
    Info<< "size after read:" << list1.size() << "\n";

    list1.set(45, list1.max_value);
    Info<< "size after write:" << list1.size() << "\n";
    Info<< "list[45]:" << list1[45] << "\n";
    list1.set(49, list1.get(100));
    report(list1);


    Info<< "\ntest copy constructor + push_back\n";
    PackedList<3> list2(list1);
    list2.push_back(4);
    Info<< "source list:\n";
    report(list1);

    Info<< "destination list:\n";
    report(list2);

    Info<< "\ntest pattern that fills all bits\n";
    PackedList<4> list3(8, 8);

    label pos = list3.size() - 1;

    list3[pos--] = list3.max_value;
    list3[pos--] = 0;
    list3[pos--] = list3.max_value;
    report(list3);

    Info<< "removed final value: " << list3.remove() << endl;
    report(list3);

    Info<<"list: " << list3 << endl;


    List<bool> list4(16, false);
    {
        // fill with some values
        forAll(list4, i)
        {
            list4[i] = i % 3;
        }

        const UList<bool>& constLst = list4;
        Info<< "\ntest operator[] const with out-of-range index\n";
        Info<< constLst << endl;
        if (constLst[100])
        {
            Info<< "[100] is true (unexpected)\n";
        }
        else
        {
            Info<< "[100] is false (expected) "
                << "list size should be unchanged (const)\n";
        }
        Info<< constLst << endl;
    }


    bitSet listb(list4);

    Info<< "copied from bool list " << endl;
    // report(listb);

    {
        labelList indices = listb.toc();

        Info<< "indices: " << indices << endl;
    }


    Info<< nl
        << "resizing: " << nl;
    {
        PackedList<1> list1(81, 1);
        PackedList<3> list3(27, 5); // ie, 101

        Info<< "initial" << nl; report(list1, true);
        Info<< "initial" << nl; report(list3, true);

        list1.resize(118, 1);
        list3.resize(37, 3);
        Info<< "extend with val" << nl; report(list1, true);
        Info<< "extend with val" << nl; report(list3, true);

        list1.resize(90, 0);
        list3.resize(30, 4);

        Info<< "contract with val" << nl; report(list1, true);
        Info<< "contract with val" << nl; report(list3, true);
    }

    {
        bitSet bits(45);

        Info<< "bits" << nl; report(bits, true);

        bits = true;
        Info<< "bits" << nl; report(bits, true);

        bits.unset(35);
        Info<< "bits" << nl; report(bits, true);

        bits.resize(39);
        Info<< "bits" << nl; report(bits, true);

        Info<< "values:" << flatOutput(bits.values()) << nl;
        Info<< "used:" << flatOutput(bits.toc()) << nl;

        bits.unset(labelRange(-15, 8));
        Info<< "bits" << nl; report(bits, true);

        bits.set(labelRange(-15, 100));
        Info<< "bits" << nl; report(bits, true);

        bits.set(labelRange(-15, 100));
        Info<< "bits" << nl; report(bits, true);

        bits.set(labelRange(150, 15));
        Info<< "bits" << nl; report(bits, true);
        bits.set(labelRange(5));
        Info<< "bits" << nl; report(bits, true);


        bits.reset();
        bits.resize(50);

        Info<< "bits" << nl; report(bits, true);
        bits.set(labelRange(4, 8));
        Info<< "bits" << nl; report(bits, true);

        bits.set(labelRange(30, 35));
        Info<< "bits" << nl; report(bits, true);

        bits.set(labelRange(80, 12));
        Info<< "bits" << nl; report(bits, true);

        bits.unset(labelRange(35, 16));
        Info<< "bits" << nl; report(bits, true);

        bits.unset(labelRange(50));
        bits.resize(100000);
        bits.set(labelRange(30, 6));
        bits[33] = false;

        Info<<"used: " << flatOutput(bits.toc()) << endl;

        Info<<"first: " << bits.find_first() << endl;
        Info<<"next: " << bits.find_next(29) << endl;
        Info<<"next: " << bits.find_next(30) << endl;
        Info<<"next: " << bits.find_next(31) << endl;

        Info<<"next: " << bits.find_next(31) << endl;

        bits.set(labelRange(80, 10));
        bits.resize(100);
        Info<< "bits" << nl; report(bits, true);

        bits.set(labelRange(125, 10));

        Info<<"next: " << bits.find_next(64) << endl;
        Info<<"used: " << flatOutput(bits.toc()) << endl;

        for (const auto pos : bits)
        {
            Info<<"have: " << pos << nl;
        }

    }

    Info<< "\n\nDone.\n";

    return 0;
}


// ************************************************************************* //
