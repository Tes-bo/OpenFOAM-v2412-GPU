/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// #define Foam_autoPtr_deprecate_setMethod

#include "autoPtr.H"
#include "labelList.H"
#include "ListOps.H"
#include "IOstreams.H"
#include "Switch.H"

#include "C7H16.H"

using namespace Foam;


// An example of bad use, since our autoPtr is too generous when being passed
// around
void testTransfer1(autoPtr<labelList> ap)
{
    // Passed in copy, so automatically removes content
    // Transfer would be nice, but not actually needed

    Info<< "recv " << Switch::name(bool(ap)) << nl;
}


// An example of good use. We are allowed to manage the memory (or not)
// and not automatically start losing things.
void testTransfer2(autoPtr<labelList>&& ap)
{
    // As rvalue, so this time we actually get to manage content
    Info<< "recv " << Switch::name(bool(ap)) << nl;
}


// Constructor from literal nullptr is implicit
template<class T>
autoPtr<T> testNullReturn1()
{
    return nullptr;
}


// Constructor from raw pointer is explicit
template<class T>
autoPtr<T> testNullReturn2()
{
    T* p = new T;

    return p;
}


template<class T>
struct DerivedList : public List<T>
{
    // Inherit constructors
    using List<T>::List;
};


template<class T>
void printInfo(const autoPtr<T>& item, const bool verbose = false)
{
    Info<< "autoPtr good:" << Switch::name(item.good())
        << " addr: " << Foam::name(item.get());

    if (verbose && item)
    {
        Info<< " content: " << item();
    }
    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    {
        auto list = autoPtr<labelList>::New(10, label(-1));

        Info<<"create: " << *list << nl;

        const labelList* plist = list;

        Info<<"pointer: " << name(plist) << nl
            <<"content: " << *plist << nl;

        Info<<"create: " << autoPtr<labelList>::New(10, label(-1))()
            << nl << nl;

        // Transfer to unique_ptr
        std::unique_ptr<labelList> list2(list.release());

        Info<<"move to unique_ptr: " << *list2 << nl;
        Info<<"old is " << Switch(bool(list)) << nl;

        autoPtr<labelList> list3(list2.release());

        Info<<"move unique to autoPtr: " << *list3 << nl;
        Info<<"old is " << Switch(bool(list2)) << nl;

        Info<< "before emplace: ";
        printInfo(list, true);

        list.emplace(4, label(-2));
        Info<< "after emplace: ";
        printInfo(list, true);

        list.emplace(2, label(-4));
        Info<< "after emplace: ";
        printInfo(list, true);
    }

    // Confirm that forwarding with move construct actually works as expected
    {
        auto source = identity(8);
        Info<<"move construct from "
            << flatOutput(source) << " @ " << name(source.cdata())
            << nl << nl;

        auto list = autoPtr<labelList>::New(std::move(source));

        Info<<"created: "
            << flatOutput(*list) << " @ " << name(list->cdata())
            << nl << nl;

        Info<<"orig: "
            << flatOutput(source) << " @ " << name(source.cdata())
            << nl << nl;
    }

    // Explicit construct Base from Derived
    {
        autoPtr<liquidProperties> liqProp
        (
            autoPtr<C7H16>::New()
        );

        Info<<"liq 1: " << liqProp() << nl << nl;
    }

    // Construct Base from Derived
    {
        autoPtr<liquidProperties> liqProp =
            autoPtr<liquidProperties>::NewFrom<C7H16>();

        Info<<"liq 2: " << liqProp() << nl << nl;
    }

    // Construct Base from Derived
    {
        const autoPtr<liquidProperties> liqProp(autoPtr<C7H16>::New());

        Info<<"liq: " << liqProp() << nl << nl;
        Info<<"liq-type: " << liqProp->type() << nl << nl;
        Info<<"type: " << typeid(liqProp.get()).name() << nl;
    }

    // Memory transfer
    {
        Info<< nl << nl;

        auto list = autoPtr<labelList>::New(identity(8));
        Info<<"forward to function from "
            << flatOutput(*list) << " @ " << name(list->cdata())
            << nl << nl;

        testTransfer2(std::move(list));

        Info<<"now have valid=" << Switch::name(bool(list));

        if (list)
        {
            Info<< nl
                << flatOutput(*list) << " @ " << name(list->cdata())
                << nl;
        }
        else
        {
            Info<< nl;
        }

        // These should fail to compile
        #if 0
        label val0 = 0;

        if (true)
        {
            val0 = list;
        }

        label val1 = 10;

        if (val1 == list)
        {
        }
        #endif
    }

    // Memory transfer
    {
        Info<< nl << nl;

        testTransfer2(autoPtr<labelList>::New(identity(8)));
    }

    // Memory transfer
    {
        Info<< nl << nl;

        auto list = autoPtr<labelList>::New(identity(8));
        Info<<"forward to function from "
            << flatOutput(*list) << " @ " << name(list->cdata())
            << nl << nl;

        testTransfer2(std::move(list));

        Info<<"now have valid=" << Switch::name(bool(list));

        if (list)
        {
            Info<< nl
                << flatOutput(*list) << " @ " << name(list->cdata())
                << nl;
        }
        else
        {
            Info<< nl;
        }
    }


    // Memory transfer
    {
        auto ptr1 = autoPtr<labelList>::New();
        auto ptr2 = autoPtr<labelList>::New();

        Info<<"ptr valid: " << bool(ptr1) << nl;

        // Refuses to compile (good!):   ptr1 = new labelList(10);

        // Does compile (good!):   ptr1 = nullptr;

        ptr1.reset(std::move(ptr2));

        autoPtr<labelList> ptr3;

        // set() method - deprecated warning?
        ptr3.set(ptr1.release());
    }


    {
        // Good this work:
        autoPtr<labelList> ptr1 = testNullReturn1<labelList>();

        // Good this does not compile:
        // autoPtr<labelList> ptr2 = testNullReturn2<labelList>();
    }

    {
        auto input1 = autoPtr<DerivedList<label>>::New(label(10), 1);
        auto input2 = autoPtr<DerivedList<scalar>>::New(label(10), 1.0);

        autoPtr<labelList> ptr1(std::move(input1));

        // Does not compile: ptr1 = std::move(input2);
        // Does not compile: ptr1 = autoPtr<List<scalar>>::New(label(10), 2);
    }


    return 0;
}


// ************************************************************************* //
