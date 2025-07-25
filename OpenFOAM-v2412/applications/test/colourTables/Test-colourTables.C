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

\*---------------------------------------------------------------------------*/

#include "colourTable.H"
#include "IOstreams.H"

using namespace Foam;

void dumpTable(const colourTable& tbl, const label n=128)
{
    Info<< tbl.table(n) << nl;
}


void dumpTable(const colourTable* tbl, const label n=128)
{
    if (tbl)
    {
        Info<< tbl->table(n) << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    // dumpTable(colourTable::ptr(colourTable::RAINBOW));
    dumpTable(colourTable::ptr(colourTable::COOL_WARM));

//     forAllConstIters(colourTable::tables(), iter)
//     {
//         Info<< nl << iter.key() << nl;
//         dumpTable(iter.val());
//     }

    Info<< "\nDone\n";

    return 0;
}


// ************************************************************************* //
