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
    Test-scalarOps

Description
    Test scalar-only ops

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "labelList.H"
#include "scalarList.H"
#include "FlatOutput.H"
#include "Tuple2.H"
#include "ops.H"
#include "scalarOps.H"
#include "vector.H"
#include "Tuple2.H"

using namespace Foam;


template<class T>
void testDivide(const List<Tuple2<T, scalar>>& list)
{
    const scalarDivideOp<T> bop;

    for (const auto& pair : list)
    {
        Info<< "num=" << pair.first()
            << " den=" << pair.second() << flush;

        Info<< " = " << bop(pair.first(), pair.second())
            << endl;
    }

    Info<< "----" << nl;

    for (const auto& pair : list)
    {
        Info<< "num=" << pair.first()
            << " den=" << pair.second() << flush;

        Info<< " = " << (pair.first() / pair.second()) << endl;
    }
}


void testModulo(const List<Tuple2<scalar, scalar>>& list)
{
    const scalarModuloOp<scalar> bop;

    for (const auto& pair : list)
    {
        Info<< "num=" << pair.first()
            << " den=" << pair.second() << flush;

        Info<< " = " << bop(pair.first(), pair.second())
            << endl;
    }

    Info<< "----" << nl;

    for (const auto& pair : list)
    {
        Info<< "num=" << pair.first()
            << " den=" << pair.second() << flush;

        Info<< " = " << std::fmod(pair.first(), pair.second()) << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    List<Tuple2<scalar, scalar>> scalars
    ({
        {10.0, 15},
        {5.0, 15},
        {5.0, 0},
    });

    List<Tuple2<vector, scalar>> vectors
    ({
        { {1,2,3}, 15},
        { {4,5,6}, 15},
        { {7,8,9}, 0},
    });


    Info<< nl << "Test scalar/scalar division" << nl;
    testDivide<scalar>(scalars);

    Info<< nl << "Test scalar/scalar modulo" << nl;
    testModulo(scalars);


    Info<< nl << "Test vector/scalar division" << nl;
    testDivide<vector>(vectors);


    //(void) compareOp<vector>()(vector::zero, vector::one);


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
