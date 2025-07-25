/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Test-boolList

Description
    Test specialized boolList functionality

\*---------------------------------------------------------------------------*/

#include "uLabel.H"
#include "boolList.H"
#include "bitSet.H"
#include "BitOps.H"
#include "FlatOutput.H"
#include "bitSetOrBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

labelList sortedToc(const UList<bool>& bools)
{
    label count = 0;
    for (const bool val : bools)
    {
        if (val) ++count;
    }

    labelList indices(count);
    count = 0;

    forAll(bools, i)
    {
        if (bools[i])
        {
            indices[count++] = i;
        }
    }
    indices.resize(count);

    return indices;
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
    const UList<bool>& bitset,
    bool showBits = false,
    bool debugOutput = false
)
{
    info(bitset);

    Info<< "values: " << flatOutput(sortedToc(bitset)) << nl;

    return Info;
}


inline Ostream& report(const UList<bool>& bools)
{
    info(bools);
    return Info;
}


inline bool compare
(
    const UList<bool>& bitset,
    const std::string& expected
)
{
    std::string has;
    has.reserve(bitset.size());

    forAll(bitset, i)
    {
        has += (bitset[i] ? '1' : '.');
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
    boolList list1(22, false);
    // Set every third one on
    forAll(list1, i)
    {
        list1[i] = !(i % 3);
    }
    Info<< "\nalternating bit pattern\n";
    compare(list1, "1..1..1..1..1..1..1..1");

    BitOps::unset(list1, labelRange(13, 20));  // In range
    Info<< "\nafter clear [13,..]\n";
    compare(list1, "1..1..1..1..1.........");

    BitOps::unset(list1, labelRange(40, 20));  // out of range
    Info<< "\nafter clear [40,..]\n";
    compare(list1, "1..1..1..1..1.........");

    BitOps::set(list1, labelRange(13, 5));  // In range
    Info<< "\nafter set [13,5]\n";
    compare(list1, "1..1..1..1..111111....");

    {
        boolList list2(5, true);
        list2.unset(2);
        list2.unset(3);

        Info<< "Test wrapper idea" << nl;

        bitSetOrBoolList wrapper(list2);

        // Use predicate, or test() method

        if (list2.test(1))
        {
            Info<< "1 is on (original)" << nl;
        }
        if (!list2.test(2))
        {
            Info<< "2 is off (original)" << nl;
        }
        if (!list2.test(100))
        {
            Info<< "100 is off (original)" << nl;
        }

        if (wrapper.test(1))
        {
            Info<< "1 is on (wrapped)" << nl;
        }
        if (!wrapper.test(2))
        {
            Info<< "2 is off (wrapped)" << nl;
        }
    }

    Info<< "\nDone" << nl << endl;
    return 0;
}


// ************************************************************************* //
