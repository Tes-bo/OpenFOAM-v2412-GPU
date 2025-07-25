/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Test the sizeof for basic types.
    Can be compiled and run without any OpenFOAM libraries.

        g++ -std=c++11 -oTest-machine-sizes Test-machine-sizes.cpp

\*---------------------------------------------------------------------------*/

#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <typeinfo>

template<class T>
void print(const char* name, bool showLimits = true)
{
    std::cout
        << "name=\"" << name << "\" sizeof=" << sizeof(T);

    if (showLimits)
    {
        std::cout
            << " \"max\"=" << std::numeric_limits<T>::max();
    }

    std::cout<< '\n';
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    std::cout<< "c++ = " << __cplusplus << '\n';
    std::cout<< "machine sizes\n---\n\n";

    print<short>("short");
    print<int>("int");
    print<long>("long");
    print<unsigned long>("unsigned long");
    print<std::size_t>("std::size_t");
    print<long long>("long long");

    print<float>("float");
    print<double>("double");
    print<long double>("long double");

    std::cout << "\n---\nEnd\n\n";

    return 0;
}


// ************************************************************************* //
