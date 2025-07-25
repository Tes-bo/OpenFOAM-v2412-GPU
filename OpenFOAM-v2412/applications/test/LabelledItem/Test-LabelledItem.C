/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Test LabelledItem (formerly 'Keyed', but that was never used)

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "edge.H"
#include "LabelledItem.H"
#include "List.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    typedef LabelledItem<edge> labelledEdge;

    List<labelledEdge> edges(10);

    forAll(edges, edgei)
    {
        auto& e = edges[edgei];

        e.insert(20-edgei);
        e.insert(edgei);

        if (!(edgei % 3))
        {
            e.setIndex(edgei);
        }
    }

    Info<< "edges: " << edges << nl;

    Foam::sort(edges);

    Info<< "sorted: " << edges << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
