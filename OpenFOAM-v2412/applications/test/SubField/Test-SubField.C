/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Test-SubField

Description
    Simple tests on SubList, SubField

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"

#include "scalarField.H"
#include "SubField.H"
#include "labelRange.H"
#include "ListOps.H"
#include <numeric>

using namespace Foam;

template<class T>
void print(const UList<T>& list)
{
    Info<< flatOutput(list) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();

    {
        List<label> ident(25);
        Foam::identity(ident, 0);

        print(ident);

        SubList<label>(ident, 10) = -10;
        print(ident);

        SubField<label>(ident, 10) = 10;
        print(ident);

        SubField<label>(ident, 10) += 10;
        print(ident);

        SubField<label>{ident, 10, 10} *= 5;
        print(ident);


        // NOTE: Need {} instead of ()
        // SubList<label>(ident) = 100;

        // GCC
        // error: conflicting declaration 'Foam::SubList<double> ident'

        // CLANG
        // warning: parentheses were disambiguated as redundant parentheses
        // around declaration of variable named 'ident' [-Wvexing-parse]

        SubList<label>{ident} = 100;
        print(ident);

        SubList<label> sub(ident);
        sub = 1;
        print(sub);

        sub.reset(ident, labelRange(4, 5)) = 5;
        print(sub);
        print(ident);

        sub.reset(ident, labelRange(14, 5)) = 15;
        print(sub);
        print(ident);

        // Cryptic, probably not a great idea to write this
        sub.reset(ident, {20, 3}) = -1;
        print(sub);
        print(ident);

        // This is also possible since we hold a concrete pointer/size
        // and not an intermediate
        Foam::identity(sub.reset(ident, 8, 8));
        print(sub);
        print(ident);
    }

    Info << "\nEnd\n";
    return 0;
}

// ************************************************************************* //
