/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Test HashTable resizing

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "HashSet.H"
#include "HashTable.H"
#include "Map.H"
#include "cpuTime.H"
#include "memInfo.H"

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

// #undef ORDERED
// #define ORDERED

using namespace Foam;

template<class T>
Ostream& printInfo(Ostream& os, const HashTable<T, T, Hash<T>>& ht)
{
    os << " (size " << ht.size() << " capacity " << ht.capacity() << ") ";
    return os;
}


template<class K, class V>
inline void insertElem
(
    #ifdef ORDERED
    std::set<K>& container,
    #else
    std::unordered_set<K, Hash<K>>& container,
    #endif
    K k,
    V v
)
{
    container.insert(k);
}


template<class K, class V>
inline void insertElem
(
    #ifdef ORDERED
    std::map<K, V>& container,
    #else
    std::unordered_map<K, V, Hash<K>>& container,
    #endif
    K k,
    V v
)
{
    container.insert(std::make_pair(k, v));
}


template<class K, class V>
inline void insertElem
(
    HashSet<K, Hash<K>>& container,
    K k,
    V v
)
{
    container.insert(k);
}


template<class K, class V>
inline void insertElem
(
    HashTable<K, V, Hash<K>>& container,
    K k,
    V v
)
{
    container.insert(k, v);
}


template<class Container>
inline void loopInsert(Container& container, const label n)
{
    for (label i = 0; i < n; i++)
    {
        insertElem(container, i, i);
    }
}


template<class Container>
inline unsigned long loopFind(const Container& container, const label n)
{
    const auto endIter = container.end();
    unsigned long sum = 0;

    for (label i = 0; i < n; i++)
    {
        if (container.find(i) != endIter)
        {
            ++sum;
        }
    }

    return sum;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    const label nLoops = 200;
    const label nFind  = 10;
    const label nElem  = 1000000;

    argList::noBanner();
    argList::addBoolOption("find", "test find");
    argList::addBoolOption("set", "test HashSet");
    argList::addBoolOption("std", "std::unordered_map or std::unordered_set");

    argList args(argc, argv);

    const bool optFnd = args.found("find");
    const bool optSet = args.found("set");
    const bool optStd = args.found("std");


    cpuTime timer;

    Info<< "insert " << nElem << " (int) elements";
    if (optFnd)
    {
        Info<< ", then find " << (nFind*nLoops) << " times\n";
    }
    else
    {
        Info<< " repeated " << nLoops << " times " << endl;
    }

    if (false)
    {
        // Verify that resizing around (0) doesn't fail
        HashTable<label, label, Hash<label>> map(32);
        printInfo(Info, map) << endl;

        map.insert(10, 1000);

        map.resize(0);
        printInfo(Info, map) << endl;

        map.resize(10);
        printInfo(Info, map) << endl;

        map.clear();
        printInfo(Info, map) << endl;

        map.resize(0);
        printInfo(Info, map) << endl;
        return 0;
    }


    if (optStd)
    {
        if (optFnd)
        {
            #ifdef ORDERED
            Info<< "using stl::set" << endl;
            std::set<label> map;
            #else
            Info<< "using stl::unordered_set" << endl;
            std::unordered_set<label, Hash<label>> map(32);
            #endif

            loopInsert(map, nElem);
            timer.resetCpuTimeIncrement();

            unsigned long sum = 0;
            for (label loopi = 0; loopi < nFind*nLoops; ++loopi)
            {
                sum += loopFind(map, nElem);
            }

            // check result (suppress compiler optimizations?)
            if (sum == 0)
            {
                Info<<"sum=0\n";
            }
        }
        else if (optSet)
        {
            #ifdef ORDERED
            Info<< "using stl::set" << endl;
            #else
            Info<< "using stl::unordered_set" << endl;
            #endif

            for (label loopi = 0; loopi < nLoops; ++loopi)
            {
                #ifdef ORDERED
                std::set<label> map;
                #else
                std::unordered_set<label, Hash<label>> map(32);
                #endif

                loopInsert(map, nElem);
            }
        }
        else
        {
            #ifdef ORDERED
            Info<< "using stl::map" << endl;
            #else
            Info<< "using stl::unordered_map" << endl;
            #endif

            for (label loopi = 0; loopi < nLoops; ++loopi)
            {
                #ifdef ORDERED
                std::map<label, label> map;
                #else
                std::unordered_map<label, label, Hash<label>> map(32);
                #endif

                loopInsert(map, nElem);
            }
        }
    }
    else
    {
        if (optFnd)
        {
            Info<< "using HashSet" << endl;

            HashSet<label, Hash<label>> map(32);

            loopInsert(map, nElem);
            timer.resetCpuTimeIncrement();

            unsigned long sum = 0;
            for (label loopi = 0; loopi < nFind*nLoops; ++loopi)
            {
                sum += loopFind(map, nElem);
            }

            // check result (suppress compiler optimizations?)
            if (sum == 0)
            {
                Info<<"sum=0\n";
            }
        }
        else if (optSet)
        {
            Info<< "using HashSet" << endl;
            for (label loopi = 0; loopi < nLoops; ++loopi)
            {
                HashSet<label, Hash<label>> map(32);

                loopInsert(map, nElem);
            }
        }
        else
        {
            Info<< "using HashTable" << endl;
            for (label loopi = 0; loopi < nLoops; ++loopi)
            {
                HashTable<label, label, Hash<label>> map(32);

                loopInsert(map, nElem);
            }
        }
    }

    Info<< timer.cpuTimeIncrement() << " s\n";
    Foam::memInfo{}.writeEntry("mem-info", Info);
    Info<< endl;

    return 0;
}

// ************************************************************************* //
