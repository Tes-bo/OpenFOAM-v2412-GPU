/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Testing of NamedEnum.
    The class is deprecated, but we may still need to support it for things
    like swak4Foam etc.

\*---------------------------------------------------------------------------*/

#include "NamedEnum.H"
#include "FlatOutput.H"

// Expect compile-time warnings

using namespace Foam;

struct testing
{
    enum class option { A, B, C, D };

    static const Foam::NamedEnum<option, 4> option1Names;
};


// All names - we have no choice with NamedEnum
template<>
const char* Foam::NamedEnum<testing::option, 4>::names[] =
{
    "a",
    "b",
    "c",
    "d"
};


const Foam::NamedEnum<testing::option, 4> testing::option1Names;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Enum 1" << nl
        << "    names:  " << testing::option1Names << nl
        << "    values: " << flatOutput(testing::option1Names.values())
        << nl << nl;

    Info<< nl
        << int(testing::option1Names["a"]) << nl
        << testing::option1Names[testing::option::A] << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
