/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test-barycentric

Description
    Some simple tests for barycentric coordinates and transforms

\*---------------------------------------------------------------------------*/

#include "barycentricTensor.H"
#include "tetrahedron.H"
#include "vectorField.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    // Tets to test
    tetPoints tetA
    (
        point(0, 0, 0),
        point(1, 0, 0),
        point(1, 1, 0),
        point(1, 1, 1)
    );

    const barycentricTensor baryT(tetA[0], tetA[1], tetA[2], tetA[3]);

    Info<< nl << "Tet: " << tetA << nl;
    Info<< "tens:" << baryT << nl;

    for
    (
        const barycentric& bary :
        List<barycentric>
        ({
            {0.25, 0.25, 0.25, 0.25},
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1},
            {0, 0, 0, 0}  // Not really valid
        })
    )
    {
        vector v(tetA.tet().barycentricToPoint(bary));
        barycentric b(tetA.tet().pointToBarycentric(v));

        Info<< nl
            << "bary: " << bary << nl
            << "vec:  " << v << nl
            // << "Vec:  " << baryT.inner(bary) << nl
            << "Vec:  " << (baryT & bary) << nl
            << "bary: " << b << nl
            // This won't work (needs a differently defined tensor)
            // << "Bary: " << (v & baryT) << nl
            ;
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
