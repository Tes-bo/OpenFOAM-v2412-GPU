/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Test-boolVector

Description
    Some simple tests for boolVector

\*---------------------------------------------------------------------------*/

#include "boolVector.H"
#include "IOstreams.H"
#include "Switch.H"

using namespace Foam;

void print(const boolVector& v)
{
    Info<< v
        << " any:" << Switch::name(v.any())
        << " all:" << Switch::name(v.all())
        << " none:" << Switch::name(v.none())
        << " on:" << v.count() << " off:" << v.count(false) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<< "boolVector" << nl
        << "  size = " << boolVector::size() << nl
        << "  contiguous = " << is_contiguous<boolVector>::value << nl
        << nl;

    {
        boolVector vec;
        Info<< "null: " << vec << nl;
    }

    Info<< "false: " << boolVector(false) << nl;
    Info<< "true: " << boolVector(true) << nl;
    Info<< "zero: " << boolVector(Zero) << nl;
    Info<< nl;

    {
        boolVector vec(1, 0, 1);
        print(vec);

        vec.flip();
        print(vec);

        vec = false;
        print(vec);

        vec = true;
        print(vec);
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
