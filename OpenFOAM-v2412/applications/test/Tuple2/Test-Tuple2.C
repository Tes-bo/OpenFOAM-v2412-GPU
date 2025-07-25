/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Application
    Test-Tuple2

Description
    Test construction, comparison etc for Tuple2 and Pair.

\*---------------------------------------------------------------------------*/

#include "labelPair.H"
#include "Tuple2.H"
#include "label.H"
#include "scalar.H"
#include "List.H"
#include "ListOps.H"
#include "ops.H"
#include "Pstream.H"
#include <functional>

using namespace Foam;

// Test for special comparison operation using compareOp
// Normal sort on label, reverse sort on scalar
struct special1
{
    typedef Tuple2<label, scalar> type;

    bool operator()(const type& a, const type& b) const
    {
        const label val = compareOp<label>()(a.first(), b.first());
        return (val == 0) ? (b.second() < a.second()) : (val < 0);
    }
};


// Test for special comparison operation using compareOp
// Normal sort on scalar, reverse sort on label
struct special2
{
    typedef Tuple2<label, scalar> type;

    bool operator()(const type& a, const type& b) const
    {
        const scalar val = compareOp<scalar>()(a.second(), b.second());
        return (val == 0) ? (b.first() < a.first()) : (val < 0);
    }
};


// Print content and info
void printTuple2(const word& f, const word& s)
{
    Info<< '(' << f << ' ' << s << ") @ "
        << name(f.data()) << ','
        << name(s.data()) << nl;
}


// Print content and info
void printTuple2(const Tuple2<word, word>& t)
{
    Info<< "tuple: " << t << " @ "
        << name(t.first().data()) << ','
        << name(t.second().data()) << nl;
}


// Print content and info
void printTuple2(const Pair<word>& t)
{
    Info<< "tuple: " << t << " @ "
        << name(t.first().data()) << ','
        << name(t.second().data()) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    typedef Tuple2<label, scalar> indexedScalar;

    Info<< "Default constructed Tuple: " << indexedScalar() << nl;
    Info<< "Default constructed Pair: "  << Pair<scalar>() << nl;

    indexedScalar t2(1, 3.2);

    Info<< "Foam::Tuple2: " << t2 << nl;

    // As list. Generated so that we have duplicate indices
    List<indexedScalar> list1(3*4);
    for (label i = 0; i < 4; ++i)
    {
        const label j = (i+1);
        const label idx = ((i % 2) ? -1 : 1) * (j);

        list1[i]   = indexedScalar(idx, (j*j));
        list1[i+4] = indexedScalar(idx, 2*j);    // duplicate index
        list1[i+8] = indexedScalar(idx+12, 2*j); // duplicate value
    }

    Info<< "Unsorted tuples:" << nl << list1 << nl;

    // Test minFirst, maxFirst functors
    {
        indexedScalar minIndexed(labelMax, Zero);
        indexedScalar maxIndexed(labelMin, Zero);

        for (const auto& item : list1)
        {
            minFirstEqOp<label>()(minIndexed, item);
            maxFirstEqOp<label>()(maxIndexed, item);
        }

        Pstream::combineReduce(minIndexed, minFirstEqOp<label>());
        Pstream::combineReduce(maxIndexed, maxFirstEqOp<label>());

        Info<< "Min indexed: " << minIndexed << nl
            << "Max indexed: " << maxIndexed << nl;
    }

    // Test minFirst, maxFirst functors
    {
        indexedScalar minIndexed(labelMax, Zero);
        indexedScalar maxIndexed(labelMin, Zero);

        for (const auto& item : list1)
        {
            minIndexed = minFirstOp<label>()(minIndexed, item);
            maxIndexed = maxFirstOp<label>()(maxIndexed, item);
        }

        Pstream::combineReduce(minIndexed, minFirstEqOp<label>());
        Pstream::combineReduce(maxIndexed, maxFirstEqOp<label>());

        Info<< "Min indexed: " << minIndexed << nl
            << "Max indexed: " << maxIndexed << nl;
    }


    Foam::sort(list1, std::less<indexedScalar>());

    Info<< "sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, std::greater<indexedScalar>());

    Info<< "reverse sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, special1());

    Info<< "special sorted tuples - sort on index, reverse on value:"
        << nl << list1 << nl;

    Foam::sort(list1, special2());

    Info<< "special sorted tuples - sort on value, reverse on index:"
        << nl << list1 << nl;


    {
        Info<< nl << nl << "Foam::Pair" << nl;

        typedef Pair<label> indexedLabel;

        indexedLabel pr(1, 3);

        Info<< "pair: "
            << pr << " => "
            << pr.first() << ' ' << pr.second() << nl;

        List<indexedLabel> list2 = ListOps::create<indexedLabel>
        (
            list1,
            [](const indexedScalar& t2)
            {
                return indexedLabel(t2.first(), t2.second());
            }
        );

        Info<< "Unsorted pairs:" << nl << list2 << nl;
    }


    {
        Info<< nl << nl << "std::pair" << nl;

        typedef std::pair<label, label> indexedLabel;

        indexedLabel pr(1, 3);

        Info<< "pair: "
            << pr << " => "
            << pr.first << ' ' << pr.second << nl;

        List<indexedLabel> list2 = ListOps::create<indexedLabel>
        (
            list1,
            [](const indexedScalar& t2)
            {
                return indexedLabel(t2.first(), t2.second());
            }
        );

        Info<< "Unsorted pairs:" << nl << list2 << nl;
    }


    {
        word word1("hello");
        word word2("word");

        Info<< "create with ";
        printTuple2(word1, word2);

        Tuple2<word> tup(std::move(word2), std::move(word1));

        printTuple2(tup);

        Info<< "input is now ";
        printTuple2(word1, word2);
    }

    {
        word word1("hello");
        word word2("word");

        Info<< "create with ";
        printTuple2(word1, word2);

        Pair<word> tup(std::move(word2), std::move(word1));

        printTuple2(tup);

        Info<< "input is now ";
        printTuple2(word1, word2);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
