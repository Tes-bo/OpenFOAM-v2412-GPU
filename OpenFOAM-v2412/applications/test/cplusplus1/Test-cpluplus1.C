/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Test miscellaneous C++ templates/functionality.

\*---------------------------------------------------------------------------*/

#include "macros.H"
#include "IOstreams.H"
#include "List.H"
#include "HashSet.H"

#include <typeinfo>
#include <type_traits>
#include <utility>

using namespace Foam;

#define PRINT_TYPEID(arg)       \
    Info<< typeid(arg).name() << " <= typeid of " << STRING_QUOTE(arg) << nl


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "various declaration types" << nl << nl;

    PRINT_TYPEID(label);
    PRINT_TYPEID(decltype(UList<label>::value_type()));
    PRINT_TYPEID(decltype(std::declval<UList<label>>().cbegin()));
    PRINT_TYPEID(decltype(*(std::declval<UList<label>>().cbegin())));
    Info<< nl;

    PRINT_TYPEID(decltype(HashTable<label>::key_type()));
    PRINT_TYPEID(decltype(HashTable<label>::value_type()));
    // Not yet: PRINT_TYPEID(decltype(HashTable<label>::mapped_type()));
    PRINT_TYPEID(decltype(std::declval<HashTable<label>>().begin()));
    PRINT_TYPEID(decltype(std::declval<const HashTable<label>>().begin()));
    PRINT_TYPEID(decltype(*(std::declval<HashTable<label>>().begin())));
    PRINT_TYPEID(decltype(*(std::declval<const HashTable<label>>().begin())));

    PRINT_TYPEID(decltype(std::declval<const HashTable<label>>().keys()));
    Info<< nl;

    PRINT_TYPEID(decltype(HashSet<label>::key_type()));
    PRINT_TYPEID(decltype(HashSet<label>::value_type()));
    // Not yet: PRINT_TYPEID(decltype(HashSet<label>::mapped_type()));
    PRINT_TYPEID(decltype(std::declval<HashSet<label>>().begin()));
    PRINT_TYPEID(decltype(std::declval<const HashSet<label>>().begin()));
    PRINT_TYPEID(decltype(*(std::declval<HashSet<label>>().begin())));
    PRINT_TYPEID(decltype(*(std::declval<const HashSet<label>>().begin())));
    Info<< nl;

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
