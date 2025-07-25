/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "argList.H"
#include "boolList.H"
#include "bitSet.H"
#include "HashSet.H"
#include "cpuTime.H"
#include <vector>
#include <unordered_set>

using namespace Foam;

#undef TEST_STD_BOOLLIST
#undef TEST_STD_UNORDERED_SET

template<class T>
void printInfo(const std::unordered_set<T, Foam::Hash<T>>& ht)
{
    Info<<"std::unordered_set elements:"
        << ht.size()
        << " buckets:" << ht.bucket_count()
        << " load_factor: " << ht.load_factor()
        << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    const label n = 1000000;
    const label nIters = 1000;

    unsigned long sum = 0;

    bitSet packed(n, true);
    boolList unpacked(n, true);

    #ifdef TEST_STD_BOOLLIST
    std::vector<bool> stdBoolList(n, true);
    #endif

    cpuTime timer;

    labelHashSet emptyHash;
    labelHashSet fullHash(1024);
    for (label i = 0; i < n; i++)
    {
        fullHash.insert(i);
    }

    Info<< "populated labelHashSet in "
        << timer.cpuTimeIncrement() << " s\n\n";

    #ifdef TEST_STD_UNORDERED_SET
    std::unordered_set<label, Foam::Hash<label>> emptyStdHash;
    std::unordered_set<label, Foam::Hash<label>> fullStdHash;
    fullStdHash.reserve(1024);
    for (label i = 0; i < n; i++)
    {
        fullStdHash.insert(i);
    }
    Info<< "populated unordered_set in "
        << timer.cpuTimeIncrement() << " s\n\n";
    #endif

    emptyHash.printInfo(Info);
    fullHash.printInfo(Info);
    #ifdef TEST_STD_UNORDERED_SET
    printInfo(emptyStdHash);
    printInfo(fullStdHash);
    #endif

    for (label iter = 0; iter < nIters; ++iter)
    {
        packed.resize(40);
        packed.shrink();
        packed.resize(n, 1);
    }
    Info<< "resize/shrink/resize:" << timer.cpuTimeIncrement() << " s\n\n";
    Info<< "packed bool size=" << packed.size() << nl;

    // Neither of these should affect the size
    packed.unset(2*n-1);
    packed.set(2*n-1, 0);
    Info<< "packed bool size=" << packed.size() << nl;

    // set every other bit on:
    Info<< "set every other bit on and count\n";
    packed.storage() = 0xAAAAAAAAu;

    // Count packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            sum += packed[i];
        }
    }

    std::cout
        << "Counting brute-force:" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;


    // Count packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        sum += packed.count();
    }

    std::cout
        << "Counting via count():" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;


    // Dummy addition
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += i + 1;
        }
    }

    std::cout
        << "Dummy loop:" << timer.cpuTimeIncrement() << " s" << nl
        << "  sum " << sum << " (sum is meaningless)" << nl;

    //
    // Read
    //

    #ifdef TEST_STD_BOOLLIST
    // Read stl
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        for (unsigned int i = 0; i < stdBoolList.size(); i++)
        {
            sum += stdBoolList[i];
        }
    }
    std::cout
        << "Reading stdBoolList:" << timer.cpuTimeIncrement() << " s" << nl
        << "  sum " << sum << nl;
    #endif

    // Read unpacked
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += unpacked[i];
        }
    }
    std::cout
        << "Reading unpacked:" << timer.cpuTimeIncrement() << " s" << nl
        << "  sum " << sum << nl;


    // Read packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            sum += packed.get(i);
        }
    }
    std::cout
        << "Reading packed using get:" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;


    // Read packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            sum += packed[i];
        }
    }
    std::cout
        << "Reading packed using reference:" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;


    // Read empty hash
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += emptyHash.found(i);
        }
    }
    std::cout
        << "Reading empty labelHashSet:" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;


    // Read full hash
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += fullHash.found(i);
        }
    }
    std::cout
        << "Reading full labelHashSet:" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;


    #ifdef TEST_STD_UNORDERED_SET
    // Read empty stl set
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += (emptyStdHash.find(i) != emptyStdHash.cend());
        }
    }
    std::cout
        << "Reading empty std::unordered_set:" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;

    // Read full stl set
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += (fullStdHash.find(i) != fullStdHash.cend());
        }
    }
    std::cout
        << "Reading full std::unordered_set:" << timer.cpuTimeIncrement()
        << " s" << nl
        << "  sum " << sum << nl;
    #endif

    Info<< "Starting write tests" << nl;

    //
    // Write
    //

    #ifdef TEST_STD_BOOLLIST
    // Read stl
    for (label iter = 0; iter < nIters; ++iter)
    {
        for (unsigned int i = 0; i < stdBoolList.size(); i++)
        {
            stdBoolList[i] = true;
        }
    }
    Info<< "Writing stdBoolList:" << timer.cpuTimeIncrement() << " s" << nl;
    #endif

    // Write unpacked
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            unpacked[i] = true;
        }
    }
    Info<< "Writing unpacked:" << timer.cpuTimeIncrement() << " s" << nl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            packed[i] = 1;
        }
    }
    Info<< "Writing packed using reference:" << timer.cpuTimeIncrement()
        << " s" << nl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            packed.set(i, 1);
        }
    }
    Info<< "Writing packed using set:" << timer.cpuTimeIncrement()
        << " s" << nl;

    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 0;
    }
    Info<< "Writing packed uniform 0:" << timer.cpuTimeIncrement()
        << " s" << nl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 1;
    }
    Info<< "Writing packed uniform 1:" << timer.cpuTimeIncrement()
        << " s" << nl;


    PackedList<3> oddPacked(n, 3);

    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 0;
    }
    Info<< "Writing packed<3> uniform 0:" << timer.cpuTimeIncrement()
        << " s" << nl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 1;
    }
    Info<< "Writing packed<3> uniform 1:" << timer.cpuTimeIncrement()
        << " s" << nl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
