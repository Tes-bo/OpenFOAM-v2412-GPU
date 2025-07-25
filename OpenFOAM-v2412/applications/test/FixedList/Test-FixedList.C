/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    Test-FixedList

Description
    Simple tests and examples for FixedList

See also
    Foam::FixedList

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "FixedList.H"
#include "Fstream.H"
#include "List.H"
#include "IPstream.H"
#include "OPstream.H"
#include <numeric>

using namespace Foam;

template<class T, unsigned N>
Ostream& printInfo(const FixedList<List<T>, N>& list)
{
    Info<< list << " addresses:";
    for (unsigned i = 0; i < N; ++i)
    {
        Info<< ' ' << name(list[i].cdata());
    }
    Info<< nl;
    return Info;
}


template<class T, unsigned N>
Ostream& printInfo
(
    const FixedList<List<T>, N>& list1,
    const FixedList<List<T>, N>& list2
)
{
    Info<< "llist1:"; printInfo(list1);
    Info<< "llist2:"; printInfo(list2);
    return Info;
}


template<class T, unsigned N>
void compileInfo()
{
    // Info<< typeid(decltype(FixedList<T, N>)).name() << nl;

    // Info<< "  holds: "
    // << typeid(decltype(FixedList<T, N>::value_type())).name() << nl;

    Info<< "max_size:"
        << FixedList<T, N>::max_size() << nl;
}


template<class FixedListType>
typename std::enable_if
<(FixedListType::max_size() == 2), bool>::type
is_pair()
{
     return true;
}


template<class FixedListType>
typename std::enable_if<(FixedListType::max_size() != 2), std::string>::type
is_pair()
{
     return "not really at all";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    argList::addBoolOption("assign");
    argList::addBoolOption("iter");
    argList::addBoolOption("swap");
    argList::addBoolOption("default", "reinstate default tests");
    argList::addBoolOption("no-wait", "test with skipping request waits");
    argList::addNote("runs default tests or specified ones only");

    #include "setRootCase.H"

    const bool optNowaiting = args.found("no-wait");

    // Run default tests, unless only specific tests are requested
    const bool defaultTests =
    (
        args.found("default")
     || args.options().empty()
     || (optNowaiting && args.options().size())
    );


    typedef FixedList<scalar,2> scalar2Type;
    typedef FixedList<label,3>  label3Type;

    // Compile-time info

    compileInfo<label, 5>();

    Info<< "pair: " << is_pair<scalar2Type>() << nl;
    Info<< "pair: " << is_pair<label3Type>() << nl;

    Info<< "max_size:" << scalar2Type::max_size() << nl;

    if (defaultTests || args.found("iter"))
    {
        Info<< nl
            << "Test iterators" << nl;

        FixedList<label, 15> ident;
        std::iota(ident.begin(), ident.end(), 0);

        // auto iter = ident.begin();
        //
        // iter += 5;
        // Info << *iter << "< " << nl;
        // iter -= 2;
        // Info << *iter << "< " << nl;

        // Don't yet bother with making reverse iterators random access
        // auto riter = ident.crbegin();

        // riter += 5;
        // Info << *riter << "< " << nl;
        // riter += 2;
        // Info << *riter << "< " << nl;

        Info<<"Ident:";
        forAllConstIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"reverse:";
        forAllReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"const reverse:";
        forAllConstReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;
    }

    if (defaultTests || args.found("swap"))
    {
        Info<< nl
            << "Test swap" << nl;

        FixedList<label, 4> list1{2, 3, 4, 5};

        Info<< "list1:" << list1
            << " hash:" << FixedList<label, 4>::hasher()(list1) << nl
            << " hash:" << Hash<FixedList<label, 4>>()(list1) << nl;

        Info<< "get<0>: " << list1.get<0>() << nl;
        Info<< "get<1>: " << list1.get<1>() << nl;
        Info<< "get<2>: " << list1.get<2>() << nl;
        Info<< "get<3>: " << list1.get<3>() << nl;
// Will not compile:  Info<< "get<4>: " << list1.get<4>() << nl;

        // Test deprecated form
        label array2[4] = {0, 1, 2, 3};
        FixedList<label, 4> list2(array2);

        Info<< "list2:" << list2
            << " hash:" << FixedList<label, 4>::hasher()(list2) << nl
            << " hash:" << Hash<FixedList<label, 4>>()(list2) << nl;


        // Using FixedList for content too
        {
            List<FixedList<label, 4>> twolists{list1, list2};
            Info<<"List of FixedList: " << flatOutput(twolists) << nl;
            sort(twolists);
            // outer-sort only
            Info<<"sorted FixedList : " << flatOutput(twolists) << nl;
        }

        Info<< "====" << nl
            << "Test swap" << nl;

        Info<< "list1: " << list1 << nl
            << "list2: " << list2 << nl;

        // Addresses don't change with swap
        Info<< "mem: "
            << name(list1.data()) << " " << name(list2.data()) << nl;

        list1.swap(list2);
        Info<< "The swap() method" << nl;
        Info<< "list1: " << list1 << nl
            << "list2: " << list2 << nl;

        Info<< "mem: "
            << name(list1.data()) << " " << name(list2.data()) << nl;

        Foam::Swap(list1, list2);
        Info<< "Foam::Swap() function" << nl;
        Info<< "list1: " << list1 << nl
            << "list2: " << list2 << nl;

        Info<< "mem: "
            << name(list1.data()) << " " << name(list2.data()) << nl;

        Info<< "====" << nl;


        Info<< nl
            << "Test of swap with other container content" << nl;

        FixedList<labelList, 4> llist1;
        FixedList<labelList, 4> llist2;

        {
            label i = 1;
            for (auto& item : llist1)
            {
                item = identity(1 + 1.5*i);
                ++i;
            }
        }

        Info<< nl
            << "initial lists" << nl;
        printInfo(llist1, llist2);

        llist2.transfer(llist1);
        Info<< nl
            << "After transfer" << nl;
        printInfo(llist1, llist2);

        llist2.swap(llist1);
        Info<< nl
            << "After swap" << nl;
        printInfo(llist1, llist2);

        llist2 = llist1;
        Info<< nl
            << "After copy assignment" << nl;
        printInfo(llist1, llist2);

        llist2 = std::move(llist1);
        Info<< nl
            << "After move assignment" << nl;
        printInfo(llist1, llist2);
    }

    Info<< nl
        << "Test construct and assignment" << nl;


    List<label> list3{0, 1, 2, 3};
    FixedList<label, 4> list4(list3);
    Info<< "list3: " << list3 << nl
        << "list4: " << list4 << nl;

    list4 = {1, 20, 3, 40};
    Info<< "list4: " << list4 << nl;

    FixedList<label, 5> list5{0, 1, 2, 3, 4};
    Info<< "list5: " << list5 << nl;

    {
        const FixedList<label, 2> indices({3, 1});
        FixedList<label, 2> list4b(list4, indices);

        Info<< "subset " << list4 << " with " << indices << " -> "
            << list4b << nl;
    }

    List<FixedList<label, 2>> list6{{0, 1}, {2, 3}};
    Info<< "list6: " << list6 << nl;

    if (UPstream::parRun())
    {
        // Fixed buffer would also work, but want to test using UList
        List<labelPair> buffer;

        DynamicList<UPstream::Request> requests;

        const label numProcs = UPstream::nProcs();
        const label startOfRequests = UPstream::nRequests();


        // NOTE: also test a mix of local and global requests...
        UPstream::Request singleRequest;

        if (UPstream::master())
        {
            // Use local requests here
            requests.reserve(numProcs);

            buffer.resize(numProcs);
            buffer[0] = labelPair(0, UPstream::myProcNo());

            for (const int proci : UPstream::subProcs())
            {
                UIPstream::read
                (
                    requests.emplace_back(),
                    proci,
                    buffer.slice(proci, 1)
                );
            }

            if (requests.size() > 1)
            {
                // Or just wait for as a single request...
                singleRequest = requests.back();
                requests.pop_back();
            }

            if (requests.size() > 2)
            {
                // Peel off a few from local -> global
                // the order will not matter (is MPI_Waitall)

                UPstream::addRequest(requests.back()); requests.pop_back();
                UPstream::addRequest(requests.back()); requests.pop_back();
            }
        }
        else
        {
            buffer.resize(1);
            buffer[0] = labelPair(0, UPstream::myProcNo());

            Perr<< "Sending to master: " << buffer << endl;

            // Capture the request and transfer to the global list
            // (for testing)

            UOPstream::write
            (
                singleRequest,
                UPstream::masterNo(),
                buffer.slice(0, 1)  // OK
                /// buffer  // Also OK
            );

            // if (singleRequest.good())
            {
                UPstream::addRequest(singleRequest);
            }
        }

        Pout<< "Pending requests [" << numProcs << " procs] global="
            << (UPstream::nRequests() - startOfRequests)
            << " local=" << requests.size()
            << " single=" << singleRequest.good() << nl;

        if (!optNowaiting)
        {
            UPstream::waitRequests(startOfRequests);
        }
        UPstream::waitRequests(requests);
        UPstream::waitRequest(singleRequest);

        Info<< "Gathered: " << buffer << endl;
    }

    return 0;
}


// ************************************************************************* //
