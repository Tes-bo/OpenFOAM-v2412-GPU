/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Test-ListOps2

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "FixedList.H"
#include "DynamicList.H"
#include "SubList.H"
#include "ListOps.H"
#include "FlatOutput.H"
#include "UPtrList.H"

using namespace Foam;


// Proof-of-concept for sorted HashTable output
// .. but yet not really convincing


// Forward declarations
template<class T, class Key, class Hash> class HashSorter;

template<class T, class Key, class Hash>
Ostream& operator<<
(
    Ostream& os,
    const HashSorter<T, Key, Hash>& sorter
);


template<class T, class Key, class Hash=Foam::Hash<Key>>
class HashSorter
{
    const HashTable<T, Key, Hash>& table;

public:

    HashSorter(const HashTable<T, Key, Hash>& ht)
    :
        table(ht)
    {}

    friend Ostream& operator<<
    (
        Ostream& os,
        const HashSorter<T, Key, Hash>& sorter
    )
    {
        const auto& tbl = sorter.table;
        const label len = tbl.size();

        // Should actually be able to get the flat entries or iterators
        // and sort that instead.

        UPtrList<const Key> keys(len);
        label count = 0;

        for (auto iter = tbl.cbegin(); iter != tbl.cend(); ++iter)
        {
            keys.set(count, &(iter.key()));
            ++count;
        }

        labelList order(identity(len));
        std::sort
        (
            order.begin(),
            order.end(),
            ListOps::less<UPtrList<const Key>>(keys)
        );

        // Size and start list delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (const label idx : order)
        {
            const auto& k = keys[idx];

            os << k << token::SPACE << tbl[k] << nl;
        }

        os << token::END_LIST;    // End list delimiter

        return os;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ListType>
void testFind(const ListType& list)
{
    Info<< nl << "list: " << flatOutput(list) << nl << endl;

    for (auto val : { 20, 35, 6, 13, 12, -12})
    {
        Info<< "lookup " << val
            << " found " << list.found(val)
            << " index " << list.find(val) << nl;
    }
}


template<class ListType>
void testMoving(ListType& list)
{
    Info<< nl << "list: " << flatOutput(list) << nl << endl;

    {
        const label i = 3;
        list.swapFirst(i);
        Info<<"swapFirst: " << i << " = " << flatOutput(list) << nl;
    }

    {
        const label i = 6;
        list.moveFirst(i);
        Info<<"moveFirst: " << i << " = " << flatOutput(list) << nl;
    }

    {
        const label i = 6;
        list.moveLast(i);
        Info<<"moveLast: " << i  << " = " << flatOutput(list) << nl;
    }

    {
        const label i = 8;
        list.swapLast(i);
        Info<<"swapLast: " << i << " = " << flatOutput(list) << nl;
    }
}


// Main program:

int main(int argc, char *argv[])
{
    List<label> list1(identity(15));
    shuffle(list1);

    FixedList<label, 15> list2(list1);
    inplaceReverseList(list2);

    DynamicList<label> list3(list1);
    inplaceReverseList(list3);

    testFind(list1);
    testFind(list2);
    testFind(list3);

    testMoving(list1);
    testMoving(list2);
    testMoving(list3);

    // Test remove
    {
        auto& list = list3;

        Info<< nl << "list: " << flatOutput(list) << nl << endl;

        list.pop_back();
        Info<<"remove = " << flatOutput(list) << nl;

        {
            const label i = 6;
            list.remove(i);
            Info<<"rmSwap: " << i << " = " << flatOutput(list) << nl;
        }

        {
            const label i = 3;
            list.remove(i);
            Info<<"rmSwap: " << i << " = " << flatOutput(list) << nl;
        }

        {
            const label i = 8;
            list.remove(i);
            Info<<"rmMove: " << i << " = " << flatOutput(list) << nl;
        }
    }


    // Test remapping
    {
        Info<< nl << "Test inplaceMapValue" << nl << nl;

        HashTable<label> input;
        typedef HashSorter<label, label> Mapper;
        typedef HashSorter<label, word> Sorter;

        for (label i=0; i < 10; ++i)
        {
            input.insert(word::printf("word%d", i), i);
        }

        Map<label> mapper;
        {
            // A mapping that does some, but not all values

            labelList rndList(identity(16));  // larger range
            shuffle(rndList);

            for (label i=0; i < 8; ++i)     // smaller sample
            {
                mapper.insert(rndList[i], 100*i);
            }
        }

        Info<< nl
            << "input: " << Sorter(input) << nl
            << "mapper: " << Mapper(mapper) << nl << nl;

        inplaceMapValue(mapper, input);

        Info<< nl << "output: " << Sorter(input) << nl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
