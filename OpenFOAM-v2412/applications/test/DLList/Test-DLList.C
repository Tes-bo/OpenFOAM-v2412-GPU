/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Test-DLList

Description
    Tests for doubly-linked lists

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "IOstreams.H"
#include "DLList.H"
#include "List.H"
#include "FlatOutput.H"
#include "ListOps.H"

using namespace Foam;

template<class T>
void printAddress(const UList<T>& list)
{
    Info<< "list addr: " << name(&list)
        << " data addr: " << name(list.cdata()) << nl;
}


template<class T>
void printAddresses(const DLList<List<T>>& sll)
{
    for (const auto& elem : sll)
    {
        printAddress(elem);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    DLList<scalar> myList{2.1, 3.4};
    myList = {2.1, 3.4, 4.3};

    for (int i = 0; i<10; i++)
    {
        myList.push_back(1.3*i);
    }

    myList.push_back(100.3);
    myList.push_back(500.3);

    Info<< "DLList<scalar>" << nl;
    Info<< nl << "flat-output: " << flatOutput(myList) << nl;

    Info<< nl << "range-for:" << nl;
    for (const auto& val : myList)
    {
        Info<< "  " << val << nl;
    }

    Info<< nl << "const_iterator:" << nl;

    forAllConstIters(myList, iter)
    {
        Info<< "  " << *iter << endl;
    }

    // Test bi-directional movement
    {
        const label n2 = myList.size()/2;

        Info<< nl << "test movement through " << flatOutput(myList) << nl;

        DLList<scalar>::const_iterator citer = myList.begin();

        for (label i=0; i<n2; ++i)
        {
            Info<< "  forward " << i << " " << *citer << nl;
            ++citer;
        }

        for (label i=0; i<n2; ++i)
        {
            Info<< "  backward " << i << " " << *citer << nl;
            --citer;
        }

        // Verify - does not compile since it uses a delete method (good!)
        DLList<scalar>::iterator iter = myList.begin();
        ++iter;
        ++iter;
        ++iter;

        Info<<" now with " << *iter << nl;
        myList.remove(iter);
        Info<<" after remove " << *iter << nl;
        ++iter;
        Info<<" after incr " << *iter << nl;
        --iter;
        --iter;

        Info<<" after 2x decr " << *iter << nl;
    }

    Info<< nl << "const_reverse_iterator:" << nl;

    forAllConstReverseIters(myList, iter)
    {
        Info<< "  " << *iter << endl;
    }


    Info<< nl << "Remove elements:" << nl;

    forAllIters(myList, iter)
    {
        Info<< "  remove " << *iter;
        myList.remove(iter);

        Info<< " => " << flatOutput(myList) << nl;
    }

    myList.push_back(500.3);
    myList.push_back(200.3);
    myList.push_back(100.3);

    Info<< nl << "Testing swapUp and swapDown:" << nl;
    Info<< " => " << flatOutput(myList) << nl;

    {
        myList.swapUp(myList.DLListBase::front());
        myList.swapUp(myList.DLListBase::back());

        Info<< nl << "swapUp => " << flatOutput(myList) << nl;
    }

    {
        myList.swapDown(myList.DLListBase::front());
        myList.swapDown(myList.DLListBase::back());

        Info<< nl << "swapDown => " << flatOutput(myList) << nl;
    }


    Info<< nl << "Transfer: " << nl;
    Info<< "original: " << flatOutput(myList) << endl;

    DLList<scalar> newList;
    newList.transfer(myList);

    Info<< nl
        << "source: " << flatOutput(myList) << nl
        << "target: " << flatOutput(newList) << nl;


    Info<< nl << "Move Construct: " << nl;

    DLList<scalar> list2(std::move(newList));

    Info<< nl
        << "in : " << flatOutput(newList) << nl
        << "out: " << flatOutput(list2) << nl;

    // Move back
    Info<< nl << "Move Assignment: " << nl;

    newList = std::move(list2);

    Info<< nl
        << "in : " << flatOutput(newList) << nl
        << "out: " << flatOutput(list2) << nl;

    // Try delete data recovery
    {
        DLList<List<label>> labList;

        for (int i = 0; i<5; i++)
        {
            labList.push_back(identity(6));
        }

        Info<< nl
            << "DLList<labelList> : " << labList << nl;

        printAddresses(labList);

        auto elem = labList.removeHead();

        Info<< " removed head" << nl;
        printAddress(elem);

        elem = labList.removeHead();

        Info<< " removed head" << nl;
        printAddress(elem);

        List<label> content1 = identity(10);

        Info<< nl
            << " move push_back ";
        printAddress(content1);

        labList.push_back(std::move(content1));

        Info<< " content " << flatOutput(content1) << nl
            << " list" << labList  << nl;

        printAddresses(labList);
        // labList.push_back(content1);
    }

    Info<< nl << "Done." << endl;

    return 0;
}


// ************************************************************************* //
