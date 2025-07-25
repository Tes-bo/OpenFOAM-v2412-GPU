/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Test behaviour of UPtrList, PtrList
\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "scalar.H"
#include "IOstreams.H"
#include "PtrDynList.H"
#include "DLPtrList.H"
#include "SLPtrList.H"
#include "plane.H"
#include "DynamicList.H"
#include "PtrListOps.H"

using namespace Foam;

class Scalar
{
    scalar data_;

public:

    static bool verbose;

    constexpr Scalar() noexcept : data_(0) {}
    Scalar(scalar val) noexcept : data_(val) {}

    ~Scalar()
    {
        if (verbose) Info<< "delete Scalar: " << data_ << endl;
    }

    scalar value() const noexcept { return data_; }
    scalar& value() noexcept { return data_; }

    autoPtr<Scalar> clone() const { return autoPtr<Scalar>::New(data_); }

    friend Ostream& operator<<(Ostream& os, const Scalar& item)
    {
        os  << item.value();
        return os;
    }
};

bool Scalar::verbose = true;


// As per
//
//     template<class T>
//     Ostream& operator<<(Ostream& os, const UPtrList<T>& list)
//
// but handle nullptr

template<class T>
Ostream& printAddr
(
    Ostream& os,
    const UPtrList<T>& list
)
{
    const label len = list.size();

    // Size and start delimiter
    os  << nl << indent << len << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    for (label i=0; i < len; ++i)
    {
        os << "addr=" << Foam::name(list.get(i)) << nl;
    }

    // End delimiter
    os << decrIndent << indent << token::END_LIST << nl;
    return os;
}


// As per
//
//     template<class T>
//     Ostream& operator<<(Ostream& os, const UPtrList<T>& list)
//
// but handle nullptr

template<class T>
Ostream& print
(
    Ostream& os,
    const UPtrList<T>& list,
    const bool debug=false
)
{
    const label len = list.size();

    // Size and start delimiter
    os  << nl << indent << len << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    for (label i=0; i < len; ++i)
    {
        const T* ptr = list.get(i);

        if (ptr)
        {
            os << *ptr << nl;
        }
        else
        {
            os << "nullptr" << nl;
        }
    }

    // End delimiter
    os << decrIndent << indent << token::END_LIST << nl;
    return os;
}


template<class T, int SizeMin>
Ostream& print
(
    Ostream& os,
    const PtrDynList<T, SizeMin>& list,
    const bool debug=false
)
{
    const label len = list.size();

    // Size and start delimiter
    os  << nl << indent << len << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    for (label i=0; i < len; ++i)
    {
        const T* ptr = list.get(i);

        if (ptr)
        {
            os << *ptr << nl;
        }
        else
        {
            os << "nullptr" << nl;
        }
    }

    if (debug)
    {
        const label cap = list.capacity();

        for (label i=len; i < cap; ++i)
        {
            const T* ptr = list.get(i);

            os << "unused " << name(ptr) << nl;
        }
    }


    // End delimiter
    os << decrIndent << indent << token::END_LIST << nl;
    return os;
}


template<class T>
Ostream& print
(
    Ostream& os,
    const UList<T*>& list
)
{
    const label len = list.size();

    // Size and start delimiter
    os  << nl << indent << len << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    for (label i=0; i < len; ++i)
    {
        const T* ptr = list[i];

        if (ptr)
        {
            os << *ptr << nl;
        }
        else
        {
            os << "nullptr" << nl;
        }
    }

    // End delimiter
    os << decrIndent << indent << token::END_LIST << nl;
    return os;
}


template<class T>
Ostream& report
(
    Ostream& os,
    const UPtrList<T>& list,
    const bool debug=false
)
{
    return print(os, list, debug);
}


template<class T, int SizeMin>
Ostream& report
(
    Ostream& os,
    const PtrDynList<T,SizeMin>& list,
    const bool debug=false
)
{
    os << "capacity=" << list.capacity() << nl;
    return print(os, list, debug);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    #if 1
    {
        DLPtrList<Scalar> llist1;
        Info<< "emplace_front: " << llist1.emplace_front(100) << nl;
        Info<< "emplace_front: " << llist1.emplace_front(200) << nl;
        Info<< "emplace_front: " << llist1.emplace_front(300) << nl;
        Info<< "emplace_back:  " << llist1.emplace_back(500)  << nl;

        Info<< "DLPtrList: " << llist1 << endl;

        Scalar::verbose = false;
        llist1.clear();
        Scalar::verbose = true;
    }
    #endif

    #if 0
    {
        DLPtrList<Scalar> llist1;
        llist1.push_front(new Scalar(100));
        llist1.push_front(new Scalar(200));
        llist1.push_front(new Scalar(300));

        auto citer = llist1.begin();

        Info<< *citer << endl;
        Info<< typeid(*citer).name() << endl;

        ++citer;
        ++citer;

        --citer;

        Info<< typeid(llist1.begin()).name() << endl;

        forAllIters(llist1, it)
        {
            Info<< typeid(*it).name() << nl
                << "reversed: " << *it << endl;
        }
        for (const auto& it : llist1)
        {
            Info<< typeid(it).name() << nl
                << "for-: " << it << endl;
        }
    }
    #endif

    // Same but as SLPtrList
    #if 0
    {
        SLPtrList<Scalar> llist1;
        llist1.push_front(new Scalar(100));
        llist1.push_front(new Scalar(200));
        llist1.push_front(new Scalar(300));

        for (const auto& it : llist1)
        {
            Info<< typeid(it).name() << nl
                << "for-: " << it << endl;
        }

        PtrList<Scalar> list1b(llist1);
        Info<< list1b << endl;
    }
    #endif

    PtrList<Scalar> list1(10);

    forAll(list1, i)
    {
        list1.set(i, new Scalar(1.3*i));
    }
    {
        auto ptr = autoPtr<Scalar>::New(10);

        Info<< "add: " << Foam::name(ptr.get());
        list1.set(0, ptr);

        Info<< "ptrlist: " << Foam::name(list1.get(0)) << nl;
        Info<< "now: " << Foam::name(ptr.get()) << nl;

        ptr = autoPtr<Scalar>::New(20);

        list1.append(ptr);
        // Delete method:  list1.push_back(ptr);
        // list1.push_back(std::move(ptr));
    }


    PtrList<Scalar> list2(15);
    Info<< "Emplace set " << list2.size() << " values" << nl;
    forAll(list2, i)
    {
        list2.emplace(i, (10 + 1.3*i));
    }

    list2.release(5);
    list2.release(10);

    {
        // Memory error (with fulldebug): const label len = (list2.size()+2);
        const label len = list2.size();
        Info<< "try_emplace " << len << " values" << nl;

        for (label i = 0; i < len; ++i)
        {
            list2.try_emplace(i, (50 + 1.3*i));
        }
    }

    PtrList<Scalar> listApp;
    for (label i = 0; i < 5; ++i)
    {
        listApp.emplace_back(1.3*i);
    }

    Info<< nl
        << "list1: " << list1 << nl
        << "list2: " << list2 << nl
        << "list-appended: " << listApp << endl;


    // Release values
    {
        DynamicList<Scalar*> ptrs;

        forAll(listApp, i)
        {
            auto old = listApp.release(i);

            if (old)
            {
                ptrs.push_back(old.release());
            }
        }

        Info<< "Released pointers from";
        print(Info, listApp) << nl;

        Info<< "Into plain list of pointers";
        print(Info, ptrs) << nl;

        PtrDynList<Scalar> newlist1(ptrs);

        Info<< "Constructed from plain list of pointers";
        print(Info, ptrs) << nl;
        print(Info, newlist1) << nl;
    }


    Info<< "indirectly delete some items via set(.., nullptr) :" << endl;
    for (label i = 2; i < 5; i++)
    {
        list1.set(i, nullptr);
    }

    {
        Info<< "range-for of list (" << list1.count_nonnull() << '/'
            << list1.size() << ") non-null entries" << nl
            << "(" << nl;
        for (const auto& item : list1)
        {
            Info<< "    " << item << nl;
        }
        Info<< ")" << nl;
    }
    {
        Info<< "iterate on non-null:" << endl;
        forAllConstIters(list1, iter)
        {
            Info<< "    " << iter.key() << " : " << iter.val() << nl;
        }
    }

    Info<< "release some items:" << endl;

    for (label i = -2; i < 5; i++)
    {
        auto old = list1.release(i);

        if (!old)
        {
            Info<< i << " was already released" << nl;
        }
    }

    Info<< "list1: ";
    print(Info, list1) << nl;

    list1.resize(list1.squeezeNull());
    Info<< "squeezed null: ";
    print(Info, list1) << nl;

    Info<< "transfer list2 -> list1:" << endl;
    list1.transfer(list2);

    Info<< "list1: " << list1 << nl
        << "list2: " << list2 << endl;

    Info<< "indirectly delete some items via setSize :" << endl;
    list1.resize(4);

    Info<< "list1: " << list1 << endl;

    {
        PtrList<Scalar> list1a(list1, false);

        Info<< "Clone constructed" << endl;
        Info<< "in:  " << list1 << nl
            << "out: " << list1a << nl
            << "addresses:" << nl;
        printAddr(Info, list1);
        printAddr(Info, list1a);
        Info<<"values:" << nl;
        print(Info, list1a);

        // This should not cause problems (ie, no deletion)
        {
            auto* ptr = &(list1a.first());
            list1a.set(0, ptr);
            Info<< "values:" << nl;
            print(Info, list1a);
        }


        PtrList<Scalar> list1b(list1a, true);

        Info<< "Reuse constructed" << endl;
        Info<< "in:  " << list1a << nl
            << "out: " << list1b << nl
            << "addresses:" << nl;
        printAddr(Info, list1a);
        printAddr(Info, list1b);


        PtrList<Scalar> list1c(list1b.clone());

        Info<< "Explicit clone()" << endl;
        Info<< "in:  " << list1b << nl
            << "out: " << list1c << nl
            << "addresses:" << nl;
        printAddr(Info, list1b);
        printAddr(Info, list1c);


        PtrDynList<Scalar> dynlist1d;
        PtrDynList<Scalar, 5> dynlist1b(list1b.clone());
        PtrDynList<Scalar, 8> dynlist1c(list1b.clone());

        Info<< "append:" << nl;
        Info<< "in: " << dynlist1b << nl
            << "in: " << dynlist1c << nl
            << "addresses:" << nl;
        printAddr(Info, dynlist1b);
        printAddr(Info, dynlist1c);

        dynlist1d.push_back(std::move(dynlist1b));
        dynlist1d.push_back(std::move(dynlist1c));

        Info<< "result:" << nl;
        print(Info, dynlist1d);

        Info<< "addresses:" << nl;
        printAddr(Info, dynlist1d);

        PtrList<Scalar> list1d;

        Info<< "append:" << nl;
        Info<< "in: " << list1b << nl
            << "in: " << list1c << nl
            << "addresses:" << nl;
        printAddr(Info, list1b);
        printAddr(Info, list1c);

        list1d.push_back(std::move(list1b));
        list1d.push_back(std::move(list1c));

        Info<< "result:" << nl;
        print(Info, list1d);

        Info<< "addresses:" << nl;
        printAddr(Info, list1d);
    }


    PtrList<Scalar> list3(std::move(list1));
    Info<< "Move constructed" << endl;

    Info<< "list1: " << list1 << nl
        << "list2: " << list2 << nl
        << "list3: " << list3 << endl;


    Info<< "Move construct:" << endl;

    PtrList<Scalar> list4(std::move(list3));

    Info<< "list3: " << list3 << nl
        << "list4: " << list4 << endl;

    Info<< "Move assign:" << endl;
    list3 = std::move(list4);

    Info<< "list3: " << list3 << nl
        << "list4: " << list4 << endl;


    Info<< "UPtrList from PtrList" << nl;

    UPtrList<Scalar> ulist1(list3);
    UPtrList<Scalar> ulist1b(list3);
    UPtrList<Scalar> ulist1c(list3);

    Info<< "ulist1: " << ulist1 << nl;
    Info<< "PtrList addresses:";
    printAddr(Info, list3);
    Info<< "UPtrList addresses:";
    printAddr(Info, ulist1);
    Info<< nl;

    ulist1c.push_back(std::move(ulist1b));

    Info<< "UPtrList append/append:";
    printAddr(Info, ulist1c);
    Info<< nl;

    {
        Info<< "UPtrList(const UPtrList&)" << nl;

        const UPtrList<Scalar>& cref = ulist1;

        UPtrList<Scalar> ulist1cp(cref);

        Info<< "src addresses:";
        printAddr(Info, cref);
        Info<< "dst addresses:";
        printAddr(Info, ulist1cp);
        Info<< nl;
    }


    Info<< "Move construct:" << endl;

    UPtrList<Scalar> ulist2(std::move(ulist1));

    Info<< "ulist1: " << ulist1 << nl
        << "ulist2: " << ulist2 << nl;

    Info<< "Copy assign:" << endl;
    ulist1 = ulist2;

    Info<< "ulist1: " << ulist1 << nl
        << "ulist2: " << ulist2 << nl;

    Info<< "Move assign:" << endl;
    ulist1 = std::move(ulist2);

    Info<< "ulist1: " << ulist1 << nl
        << "ulist2: " << ulist2 << nl;

    // Test iterator random access
    #if (OPENFOAM <= 2212)
    {
        auto iter1 = ulist1.begin();
        auto iter2 = iter1 + 3;

        Info<< "begin:" << *iter1 << " (+3):" << *iter2 << nl;
        Info<< "diff= " << (iter1 - iter2) << nl;
        Info<< "iter[2]=" << iter1[2] << nl;
        Info<< "iter1 < iter2 : " << (iter1 < iter2) << nl;
        Info<< "iter1 >= iter2 : " << (iter1 >= iter2) << nl;

        Info<< "->" << iter1->value() << nl;
        Info<< "*"  << (*iter1).value() << nl;
        Info<< "()" << iter1().value() << nl;
    }
    #endif

    PtrList<plane> planes;
    planes.emplace_back(vector::one, vector::one);
    planes.emplace_back(vector(1,2,3), vector::one);

    Info<< nl << "appended values" << nl;
    for (const plane& p : planes)
    {
        Info<< "    plane " << p << endl;
    }

    Info<< "Testing PtrDynList" << nl;

    PtrDynList<plane> dynPlanes;

    {
        dynPlanes.emplace_back(vector::one, vector::one);
        dynPlanes.emplace_back(vector(1,2,3), vector::one);
        dynPlanes.push_back(nullptr);

        dynPlanes.set(6, new plane(vector(2,2,1), vector::one));
        dynPlanes.set(10, new plane(vector(4,5,6), vector::one));

        Info<< "emplaced[12]: "
            << dynPlanes.emplace(12, vector(3,2,1), vector::one) << endl;

        dynPlanes.emplace_back(Zero, vector::one);
    }

    Info<< nl << "PtrDynList: ";
    report(Info, dynPlanes, true);

    dynPlanes.resize(9);

    Info<< nl << "resize()";
    report(Info, dynPlanes, true);

    dynPlanes.clear();
    Info<< "clear()" << nl;
    report(Info, dynPlanes);

    Info<< "now append again" << endl;
    {
        dynPlanes.emplace_back(vector::one, vector::one);
        dynPlanes.emplace_back(vector(1,2,3), vector::one);

        dynPlanes.emplace(5, vector(2,2,1), vector::one);
    }

    report(Info, dynPlanes, true);

    {
        // No clone for plane - do manual copy
        PtrList<plane> stdPlanes(dynPlanes.size());

        forAll(dynPlanes, i)
        {
            const plane* pln = dynPlanes.get(i);
            if (pln)
            {
                stdPlanes.set(i, new plane(*pln));
            }
        }

        report(Info, stdPlanes);
        printAddr(Info, stdPlanes);

        stdPlanes.resize(stdPlanes.squeezeNull());

        Info<< "After pruning nullptr entries" << endl;
        printAddr(Info, stdPlanes);
    }

    dynPlanes.resize(dynPlanes.squeezeNull());

    Info<< "After pruning nullptr entries" << endl;
    report(Info, dynPlanes, true);


    {
        PtrDynList<plane> dynPlanes2;

        dynPlanes2.emplace_back(vector::one, vector::one);
        dynPlanes2.emplace_back(vector(1,2,3), vector::one);
        dynPlanes2.push_back(nullptr);

        dynPlanes2.emplace(6, vector(2,2,1), vector::one);
        dynPlanes2.emplace(10, Zero, vector::one);

        labelList order;
        sortedOrder(dynPlanes2, order);
        Info<< "sorted order: " << flatOutput(order) << nl;

        sortedOrder(dynPlanes2, order, UPtrList<plane>::greater(dynPlanes2));
        Info<< "sorted order: " << flatOutput(order) << nl;

        // Shuffle
        shuffle(dynPlanes2);
        Info<< "Shuffled" << endl;
        report(Info, dynPlanes2, false);

        // Reverse sort
        sort
        (
            dynPlanes2,
            [](const plane& a, const plane& b) { return (b < a); }
        );
        Info<< "Reverse sorted" << endl;
        report(Info, dynPlanes2, false);

        // Forward sort
        sort(dynPlanes2);
        Info<< "Sorted" << endl;
        report(Info, dynPlanes2, false);

        // Reverse pointer list - not yet available or needed!
        /// reverse(dynPlanes2);
        /// Info<< "Reversed" << endl;
        /// report(Info, dynPlanes2, false);

        // dynPlanes2.squeezeNull();

        Info<< "Append" << endl;
        report(Info, dynPlanes2, false);

        dynPlanes.push_back(std::move(dynPlanes2));

        Info<< "Result" << endl;
        report(Info, dynPlanes, false);

        Info<< "From" << endl;
        report(Info, dynPlanes2, false);
    }

    Info<< "free()" << endl;

    dynPlanes.free();
    report(Info, dynPlanes, true);

    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
