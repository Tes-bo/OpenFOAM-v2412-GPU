/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test some string_view functionality

\*---------------------------------------------------------------------------*/

#include "string.H"
#include "IOstreams.H"
#include "List.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Compiled with C++ " << __cplusplus;
    #if __cplusplus >= 201703L
    Info<< " - has std::string_view" << nl << nl;
    #else
    Info<< " - NO std::string_view" << nl << nl;
    #endif

    // basics
    {
        for
        (
            const auto& cstr
          :
            {
                "abcdef"
            }
        )
        {
            const auto len = strlen(cstr);

            Info<< nl
                << "input: <" << cstr << '>'
                << " type: " << typeid(cstr).name() << " len:" << len << nl;

            #if __cplusplus >= 201703L
            Info<< "    view: " << std::string_view(cstr) << nl;
            #endif

            Info<< "    span: "
                << stdFoam::span<const char>(cstr, len) << nl;
            Info<< "    span: "
                << stdFoam::span<char>(const_cast<char*>(cstr), len) << nl;
        }
    }

    // This should fail to compile:
    #if 0
    {
        labelList values(identity(4));

        Info<< "values: " << values << nl;

        Info<< "    span: "
            << stdFoam::span<label>(values.data(), values.size()) << nl;
    }
    #endif

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
