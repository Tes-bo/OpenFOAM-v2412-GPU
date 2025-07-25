/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Test-Map

Description

\*---------------------------------------------------------------------------*/

#include "Map.H"
#include <map>
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Map<bool> map1
    {
        {1, true}, {2, false}, {3, true}, {4, false}, {5, true}
    };

    // Taking a const iterator from find does not work!
    // Also, fails later on op==
    Map<bool>::const_iterator map1Iter = map1.cfind(5);

    // Same, but with non-const access
    // Map<bool>::iterator map1Iter = map1.find(5);

    if (!map1Iter.good())
    {
        Info<< "not found" << endl;
    }
    else
    {
        Info<< "5 is " << map1Iter.val() << endl;
    }

    // Repeat with std::map
    Info<< "Same with STL" << endl;

    std::map<label, bool> stdmap1
    {
        {1, true}, {2, false}, {3, true}, {4, false}, {5, true}
    };

    std::map<label, bool>::const_iterator stdmap1Iter = stdmap1.find(5);

    if (stdmap1Iter == stdmap1.cend())
    {
        Info<< "not found" << endl;
    }
    else
    {
        Info<< "5 is " << stdmap1Iter->second << endl;
    }


    Info<<"test move construct" << nl;
    Map<bool> map2(std::move(map1));
    Map<bool> map3;

    std::map<label, bool> stdmap2(std::move(stdmap1));
    std::map<label, bool> stdmap3;

    Info<<"map1: " << map1 << nl
        <<"map2: " << map2 << nl;

    Info
        <<"stdmap1: " << stdmap1.size() << nl
        <<"stdmap2: " << stdmap2.size() << nl;


    Info<<"test move assign" << nl;
    map3 = std::move(map2);
    stdmap3 = std::move(stdmap2);

    Info<<"map2: " << map2 << nl
        <<"map3: " << map3 << nl;

    Info
        <<"stdmap2: " << stdmap2.size() << nl
        <<"stdmap3: " << stdmap3.size() << nl;


    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
