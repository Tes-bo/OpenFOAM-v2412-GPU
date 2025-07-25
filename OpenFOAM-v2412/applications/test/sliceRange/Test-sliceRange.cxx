/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Test slice range

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelList.H"
#include "FixedList.H"
#include "sliceRange.H"
#include "IndirectList.H"
#include "IndirectSubList.H"
#include "SliceList.H"
#include "Random.H"

using namespace Foam;


typedef FixedList<label, 3> sliceCoeffs;

void printInfo(const sliceCoeffs& coeffs)
{
    sliceRange range(coeffs);

    Info<< nl
        << "coeffs: " << coeffs << nl
        << "range:  " << range << nl
        << "min:    " << range.min() << nl
        << "*begin  " << *range.begin() << nl
        << "max:    " << range.max() << nl
        << "*end    " << *range.end() << nl
        << "values: " << flatOutput(range.labels()) << nl;

    Info<< "for  :";
    for (const label val : range)
    {
        Info<< ' ' << val;
    }
    Info<< nl;

    if (range.empty())
    {
        Info<< "empty"<< nl;
        return;
    }

    Info<< "mid-point [" << (range.size()/2) << "] = "
        << range[range.size()/2] << nl;

    const auto endIter = range.cend();

    for (const label i : {label(-1), (range.size()/2), range.size()})
    {
        const auto iter = range.at(i);

        Info<< "at(" << i << ") = " << *iter;

        if (iter == endIter)
        {
            Info<< " (out-of-range)";
        }

        Info<< nl;
    }

    // Copy construct
    sliceRange range2(range);

    // Copy assign
    range2 = range;
}


void printForLoop(const sliceRange& range)
{
    Info<< "for " << range << nl
        << "  >";

    for (const label val : range)
    {
        Info<< ' ' << val;
    }

    Info<< nl;
}


template<class IteratorType>
void printIteratorTest(IteratorType& iter)
{
    const auto iter2 = (iter - 5);
    const auto iter3 = (iter + 5);

    // Info<< typeid(iter).name() << nl;

    Info<< "begin: " << *iter++;
    Info<< " next: " << *iter;
    Info<< " next: " << *(++iter);
    Info<< " [5]: " << iter[5];
    Info<< " +10: " << *(iter + 10);
    Info<< " -10: " << *(iter - 10);
    Info<< nl;

    Info<< "compare: " << *iter2 << " and " << *iter3 << nl;
    Info<< "   == " << (iter2 == iter3) << nl;
    Info<< "   <  " << (iter2 < iter3) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();

    argList args(argc, argv, false, true);


    for
    (
        sliceCoeffs coeffs :
        {
            sliceCoeffs{25, 8, 2},
            sliceCoeffs{15, 5, 3},
            sliceCoeffs{15, -5, 2},
        }
    )
    {
        printInfo(coeffs);
    }

    // Some iterator tests
    {
        sliceRange range(25, 8, 3);

        auto iter1 = range.begin();
        Info<< nl << "Forward iterator for " << range << nl;
        printIteratorTest(iter1);

        auto iter2 = range.rbegin();
        Info<< nl << "Reverse iterator for " << range << nl;
        printIteratorTest(iter2);


        // Change the start, size
        range.start() = 5;
        range.size() = 4;

        iter1 = range.begin();

        Info<< nl << "Forward iterator for " << range << nl;
        printIteratorTest(iter1);
    }


    // Generator
    {
        sliceRange range(25, 8, 3);

        Info<< nl << "Generator for " << range << nl;

        auto gen = range.generator();

        for (label i=0; i < 10; ++i)
        {
            Info<< ' ' << gen();
        }
        Info<< " ... " << nl;

        Info<< nl << "range-for:" << nl;
        for (const auto& val : range)
        {
            Info<< ' ' << val;
        }
        Info<< nl;

        Info<< nl << "reversed:" << nl;
        forAllConstReverseIters(range, iter)
        {
            Info<< ' ' << *iter;
        }
        Info<< nl;
    }

    // Sliced lists
    {
        List<scalar> list1(100);

        Random rnd(1234);

        for (scalar& val : list1)
        {
            val = 100 * rnd.sample01<scalar>();
        }

        Info<< nl << "Random list: " << flatOutput(list1) << nl;

        {
            IndirectSubList<scalar> sublist1(list1, labelRange(10));

            Info<< nl << "SubList: " << sublist1.addressing() << " = "
                << flatOutput(sublist1) << nl;

            // Adjust addressing
            sublist1.addressing().reset(5, 8);

            // This should resolve as a no-op
            sublist1 = sublist1;

            Info<< "SubList: " << sublist1.addressing() << " = "
                << flatOutput(sublist1) << nl;
        }


        SliceList<scalar> slice1(list1, sliceRange(0, 15, 3));

        Info<< nl << "slicing with: " << slice1.addressing() << nl;

        Info<< nl << "Sliced list: " << flatOutput(slice1) << nl;

        for (scalar& val : slice1)
        {
            val = -val;
        }

        // Changed list via slice:
        Info<< nl << "Changed via slice: " << flatOutput(list1) << nl;

        // Some indirect list

        IndirectList<scalar> indlist
        (
            list1,
            identity(slice1.size(), list1.size()-slice1.size())
        );

        Info<< nl << "Indirect slice: " << flatOutput(indlist) << nl;

        indlist = 1000;
        Info<< nl << "zeroed slice: " << flatOutput(indlist) << nl;

        slice1 = indlist;

        Info<< nl << "self-copy: " << flatOutput(list1) << nl;

        slice1 = 100000;

        Info<< nl << "set values: " << flatOutput(slice1) << nl
            << " = " << flatOutput(list1) << nl;
    }


    // For loops
    {
        Info<< nl << "Test for loops" << nl;

        printForLoop(sliceRange(25, 8, -2));
        printForLoop(sliceRange(10, 3, 0));
        printForLoop(sliceRange(10, 3, 2));
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
