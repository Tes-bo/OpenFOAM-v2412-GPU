/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-IntRange

Description
    Test integer range
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelPair.H"
#include "IntRange.H"
#include "SpanStream.H"

using namespace Foam;

template<class T>
void printInfo(const IntRange<T>& range)
{
    Info<< "    min  " << range.min() << nl
        << "    max  " << range.max() << nl
        << "    size   " << range.size() << nl
        << "begin end  " << *range.cbegin() << ' ' << *range.cend() << nl;

    Info<< "rbegin rend " << *range.rbegin() << ' ' << *range.rend() << nl;
}


template<class T>
void printValues(const IntRange<T>& range)
{
    Info<< range.size() << "(";
    for (const label val : range)
    {
        Info<< ' ' << val;
    }
    Info<< " )";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();

    argList args(argc, argv, false, true);

    // Fails static_assert
    /// Info<< "Default construct float: " << IntRange<float>().size() << nl;

    // Does not really make sense, but...
    /// Info<< "Default construct bool: " << IntRange<bool>().size() << nl;


    typedef IntRange<int> intRange;

    Info<< "Default construct int32_t: " << IntRange<int32_t>() << nl
        << "Default construct int64_t: " << IntRange<int64_t>() << nl;

    Info<< "  one: " << intRange(10) << nl
        << "  two: " << intRange(5, 10) << nl;

    // Read from stream
    {
        ICharStream is("(10 100)");
        intRange range;
        is >> range;

        Info<< "From stream int32_t: " << range << nl;
    }

    for
    (
        const labelPair& pr
      : labelPairList
        {
            {10, 2},
            {2, 3},
            {100, 0}
        }
    )
    {
        intRange range(pr.first(), pr.second());

        Info<< "range: " << range
            << (range ? " non-empty" : " Empty") << nl;
    }

    {
        const intRange range1(3, 16);

        auto begIter = range1.begin();
        auto endIter = range1.end();
        auto midIter = range1.at(range1.size()/2);

        Info<< "iterator tests on " << range1 << nl;
        Info<< "beg = " << *begIter << nl
            << "end = " << *endIter << nl
            << "mid = " << *midIter << nl
            << "end - beg = " << (endIter - begIter) << nl;

        Info<< "distance: " << std::distance(begIter, endIter) << nl;

        Info<< "beg + 10 = " << *(begIter + 10) << nl
            << "beg[100] = " << begIter[100] << nl;

// Info<< "10 + beg = " << *(10 + begIter) << nl;
// Will not work:
// Avoid this definition since it participates in too many resolution
// attempts and ruins everything.

        std::swap(begIter, endIter);
        Info<< "after iterator swap" << nl
            << "beg = " << *begIter << nl
            << "end = " << *endIter << nl;

        auto rbegIter = range1.rbegin();
        auto rendIter = range1.rend();

        Info<< nl
            << "reverse beg = " << *rbegIter << nl
            << "reverse end = " << *rendIter << nl
            << "reverse end - beg = " << (rendIter - rbegIter) << nl;

        Info<< "reverse beg + 10 = " << *(rbegIter + 10) << nl
            << "reverse beg[100] = " << rbegIter[100] << nl;

        std::swap(rbegIter, rendIter);
        Info<< "after iterator swap" << nl
            << "reverse beg = " << *rbegIter << nl
            << "reverse end = " << *rendIter << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
