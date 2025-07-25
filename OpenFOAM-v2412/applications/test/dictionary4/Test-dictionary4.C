/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Test-dictionary4

Description
    Test expansion
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"
#include "Pair.H"
#include "stringOps.H"

using namespace Foam;

void addToDict
(
    dictionary& dict,
    std::initializer_list<Pair<std::string>> entries
)
{
    for (const Pair<std::string>& e : entries)
    {
        dict.add(word(e.first()), string(e.second()));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    // Test expansion
    {
        dictionary dict;


        addToDict
        (
            dict,
            {
                {"fileDir", "<case>/triSurface"},
                {"fileBase", "input"},
                {"fileExt", ".stl"},
                {"fileName", "${fileDir}/${fileBase}$fileExt"}
            }
        );

        Info<< "dict" << dict << nl;

        string str("fileName = <$fileName>");

        Info<< str.c_str() << nl;

        stringOps::inplaceExpand(str, dict, true, true);

        Info<< str.c_str() << nl;
    }

    return 0;
}


// ************************************************************************* //
