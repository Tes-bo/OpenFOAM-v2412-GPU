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

Application
    Test-partition-sort

Description
    Test behaviour of std::partition, etc.

\*---------------------------------------------------------------------------*/

#include "IndirectList.H"
#include "SortList.H"
#include "Random.H"

#include <algorithm>
#include <iterator>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    labelList input(40);
    labelList order(Foam::identity(input.size()));

    UIndirectList<label> sorted(input, order);

    std::generate_n
    (
        input.begin(),
        input.size(),
        Random::uniformGeneratorOp<label>(123456, -5, 15)
    );

    Info<< "input: "; input.writeList(Info) << nl;
    //Info<< "input: "; sorted.writeList(Info) << nl;

    // Partition with 0/+ve values first, -ve values second
    auto part2 =
        std::stable_partition
        (
            order.begin(),
            order.end(),
            [&](const label a) { return (input[a] >= 0); }
        );

    Info<< "partitioned at: "
        << label(std::distance(order.begin(), part2)) << nl;


    //Info<< "partitioned: "; input.writeList(Info) << nl;
    Info<< "partitioned: "; sorted.writeList(Info) << nl;

    // Sort -ve values with increasing distance away from 0, so a reverse sort
    std::stable_sort
    (
        part2,
        order.end(),
        [&](const label a, const label b) { return (input[b] < input[a]); }
    );

    //Info<< "sorted: "; input.writeList(Info) << nl;
    Info<< "sorted: "; sorted.writeList(Info) << nl;
    Info<< "order: "; order.writeList(Info) << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
