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
    Test-HashPtrTable

Description
    Tests for HashPtrTable base functionality

\*---------------------------------------------------------------------------*/

#include <memory>
#include <iostream>
#include "autoPtr.H"
#include "refPtr.H"
#include "tmp.H"
#include "PtrMap.H"
#include "primitiveFields.H"

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

    const scalar& value() const noexcept { return data_; }
    scalar& value() noexcept { return data_; }

    friend Ostream& operator<<(Ostream& os, const Scalar& item)
    {
        os  << item.value();
        return os;
    }
};

bool Scalar::verbose = true;


template<class T>
void printTable(const HashPtrTable<T>& table)
{
    Info<< table.size() << nl << "(" << nl;

    forAllConstIters(table, iter)
    {
        const T* ptr = iter.val();
        Info<< iter.key() << " = ";
        if (ptr)
        {
            Info<< *ptr << " (" << name(ptr) << ")";
        }
        else
        {
            Info<< "nullptr";
        }
        Info<< nl;
    }

    Info<< ")" << endl;

    // Iterate across values, with for-range
    Info<< "values (";
    for (const auto& ptr : table)
    {
        Info<< ' ';
        if (ptr)
        {
            Info<< *ptr;
        }
        else
        {
            Info<< "nullptr";
        }
    }
    Info<< " )" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    HashPtrTable<double> myTable;
    myTable.set("abc", new double(42.1));
    myTable.set("def", nullptr);
    myTable.set("pi", new double(3.14159));
    myTable.set("natlog", new double(2.718282));
    myTable.insert("sqrt2", autoPtr<double>::New(1.414214));
    myTable.insert("euler", autoPtr<double>::New(0.577216));
    myTable.set("def_0", nullptr);
    myTable.emplace_set("def_1", 123);
    myTable.emplace_set("def_2", 456);

    HashTable<std::unique_ptr<double>> myTable1;

    myTable1.set("abc", std::unique_ptr<double>(new double(42.1)));
    myTable1.set("pi", std::unique_ptr<double>(new double(3.14159)));
    myTable1.set("natlog", std::unique_ptr<double>(new double(2.718282)));

    HashTable<autoPtr<double>> myTable2;

    myTable2.set("abc", autoPtr<double>(new double(42.1)));
    myTable2.set("pi", autoPtr<double>(new double(3.14159)));
    myTable2.set("natlog", autoPtr<double>(new double(2.718282)));
    myTable2.insert("sqrt2", autoPtr<double>::New(1.414214));

    Info<< "Initial table" << nl;
    printTable(myTable);

    myTable.try_emplace("def_0", 1000);  // was nullptr, now value
    myTable.try_emplace("def_1", 1001);  // no-op
    myTable.try_emplace("def_2", 1002);  // no-op;
    myTable.try_emplace("def_3", 1003);  // was non-existent, now value

    Info<< "after try_emplace" << nl;
    printTable(myTable);

    Info<< "print" << nl;
    Info<< myTable2 << nl;

    {
        auto iter2 = myTable2.find("pi");
        Info<< nl << "Got pi=";
        if (iter2.good())
        {
            Info<< **iter2 << nl;
        }
        else
        {
            Info<< "not-found" << nl;
        }
    }

    HashPtrTable<double> copy(myTable);

    // Info<< copy << endl;
    printTable(copy);
    Info<< copy << endl;

    Info<<"\nerase some existing and non-existing entries" << nl;

    auto iter = myTable.find("pi");
    myTable.erase(iter);

    iter = myTable.find("unknownKey");
    myTable.erase(iter);

    myTable.erase("abc");
    myTable.erase("unknownKey");

    myTable.release("euler");

    Info<< "After erasure" << nl;
    printTable(myTable);

    HashPtrTable<double> moved(std::move(copy));

    Info<< nl << "test movable" << nl;
    Info<<"input:" << nl;
    printTable(copy);

    Info<<"output:" << nl;
    printTable(moved);

    HashPtrTable<double> other;

    Info<<"move assign" << nl;

    other = std::move(moved);
    printTable(other);

    Info<<"old" << nl;
    printTable(moved);


    Info<< "Verifying deletion characteristics" << nl;
    {
        HashPtrTable<Scalar> tbl;
        tbl.set("abc", new Scalar(42.1));
        tbl.set("def", nullptr);
        tbl.set("pi", new Scalar(3.14159));
        tbl.set("natlog", new Scalar(2.718282));
        tbl.insert("sqrt2", autoPtr<Scalar>::New(1.414214));

        Info<< "Table: " << tbl << nl;

        Info<< nl << "Check exists, non-null" << nl;

        for (const word& k : { "abc", "foo", "pi" })
        {
            Info<< "    " << k << ' ';

            const auto* inspect = tbl.get(k);

            if (inspect)
            {
                Info<< *inspect << nl;
            }
            else
            {
                Info<< "(null)" << nl;
            }
        }

        Info<< nl << "... overwrite again" << nl;

        tbl.set("abc", new Scalar(42.1));
        tbl.set("def", nullptr);
        tbl.set("pi", new Scalar(3.14159));
        tbl.set("natlog", new Scalar(2.718282));

        tbl.emplace("other", 15.6);

        Info<< "Table: " << tbl << nl;

        Info << nl << "Test emplace and emplace_set" << nl;

        tbl.emplace("abc", 100);
        tbl.emplace_set("def", 100);
        tbl.emplace_set("other", 100);

        Info<< "Table: " << tbl << nl;
    }

    {
        PtrMap<scalarField> fields;

        {
            const label patchi = 2;

            scalarField fld1(patchi, 5.0);
            scalarField fld2(patchi, 8.0);

            // assign from tmp<>
            fields.set( patchi, (fld1 * fld2));
        }

        {
            const label patchi = 3;

            scalarField fld1(patchi, 6.0);

            // From tmp (clone)
            fields.set(patchi, tmp<scalarField>(fld1));
        }

        {
            const label patchi = 4;

            // From refPtr
            fields.set(patchi, refPtr<scalarField>::New(patchi, 10.0));
        }

        Info<< nl
            << "PtrMap:" << nl
            << fields << endl;
    }

    Info<< "\nEnd" << nl;

    return 0;
}

// ************************************************************************* //
