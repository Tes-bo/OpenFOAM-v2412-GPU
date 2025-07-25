/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test-write-wrapped-string

Description
    Simple tests for wrapped strings

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "stringOps.H"

using namespace Foam;

void print(const std::string& str, std::size_t width, std::size_t indent=0)
{
    auto& os = Info();

    os  << nl
        << "string[" << str.size() << "]" << nl
        << str.c_str() << "<<<<" << nl
        << "indent:" << indent << " width:" << width << endl;

    for (size_t i = 0; i < width; ++i)
    {
        os << '=';
    }
    os << endl;

    stringOps::writeWrapped(os, str, width, indent);

    for (size_t i = 0; i < width; ++i)
    {
        os << '=';
    }
    os << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    #include "setRootCase.H"

    {
        string test =
            "123456789-12345\n\n"
            "6789-12\t"
            "xyz3456789-1234    56789-123456789-";

        print(test, 10, 4);
    }

    {
        string test = "ABCDEFGHI";
        print(test, 10, 4);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
