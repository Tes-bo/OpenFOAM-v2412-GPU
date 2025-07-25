/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Test-nullObject

Description
    Tests of nullObject

\*---------------------------------------------------------------------------*/

#include "nullObject.H"
#include "List.H"
#include "HashSet.H"
#include "faceList.H"
#include "pointField.H"
#include "globalIndex.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

class SimpleClass
{
public:

    //- Default construct
    SimpleClass() = default;
};


template<class T>
void printInfo(const UList<T>& list)
{
    std::cout
        << nl
        << "List : addr: " << name(&list)
        << " (null: " << isNull(list) << ")" << nl
        << "    size: " << list.size() << " empty: " << list.empty() << nl
        << "    data: " << name(list.cdata())
        << " begin=" << name(list.begin())
        << " end=" << name(list.end()) << nl;

    Info<< list << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    // Test pointer and reference to a class
    auto ptrToClass = autoPtr<SimpleClass>::New();
    auto& refToClass = ptrToClass.ref();

    std::cout
        << "nullObject addr=" << name(&(nullObjectPtr)) << nl
        << "  sizeof(nullObject) = " << sizeof(NullObject::nullObject) << nl
        << "  sizeof(void*) = " << sizeof(void*) << nl
        << "  sizeof(labelList) = " << sizeof(labelList) << nl
        << "  sizeof(wordHashSet) = " << sizeof(wordHashSet) << nl << nl;

    std::cout
        << "nullObject" << nl
        << "  pointer:" << name(nullObjectPtr->pointer()) << nl
        << "  value:"   << nullObjectPtr->value() << nl << nl;

    if (notNull(ptrToClass.get()))
    {
        Info<< "Pass: ptrToClass is not null" << nl;
    }
    else
    {
        Info<< "FAIL: ptrToClass is null" << nl;
    }

    if (notNull(refToClass))
    {
        Info<< "Pass: refToClass is not null" << nl;
    }
    else
    {
        Info<< "FAIL: refToClass is null" << nl;
    }


    // Test pointer and reference to the nullObject

    const SimpleClass* ptrToNull = NullObjectPtr<SimpleClass>();
    const SimpleClass& refToNull = (*ptrToNull);

    if (isNull(ptrToNull))
    {
        Info<< "Pass: ptrToNull is null" << nl;
    }
    else
    {
        Info<< "FAIL: ptrToNull is not null" << nl;
    }

    if (isNull(refToNull))
    {
        Info<< "Pass: refToNull is null" << nl;
    }
    else
    {
        Info<< "FAIL: refToNull is not null" << nl;
    }


    // Test List casting
    {
        labelList list1;
        labelList list2({1, 2, 3});

        printInfo(list1);
        printInfo(list2);
        printInfo(labelList::null());

        printInfo(faceList::null());
        printInfo(pointField::null());
    }

    // Test swallow assigment (like std::ignore)
    // Looks pretty ugly though!

    NullObject::nullObject = "hello world";
    NullObject::nullObject = Foam::identity(5);


    {
        const auto& gi = globalIndex::null();
        Info<< "globalIndex::null() => "
            << " empty: " << gi.empty()
            << " nProcs: " << gi.nProcs()
            << " total-size: " << gi.totalSize() << nl;

        // Even this works
        Info<< "   offsets: " << gi.offsets() << nl;
    }

    Info<< nl;

    return 0;
}

// ************************************************************************* //
