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

Description
    Tests for tensor2D and vector2D

\*---------------------------------------------------------------------------*/

#include "tensor2D.H"
#include "vector2DField.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    vector2D v1(1, 2), v2(3, 4);
    tensor2D t3 = v1*v2;

    Info<< v1 << "*" << v2 << " = " << t3 << endl;

    {
        Info<< "rows:" << nl;
        for (direction i=0; i < 2; ++i)
        {
            Info<< "  (" << i << ") = " << t3.row(i) << nl;
        }
    }

    {
        Info<< "cols:" << nl;
        for (direction i=0; i < 2; ++i)
        {
            Info<< "  (" << i << ") = " << t3.col(i) << nl;
        }
        Info<< "col<0> = " << t3.col<0>() << nl;
        Info<< "col<1> = " << t3.col<1>() << nl;
        // Compilation error:  Info << "col<3> = " << t3.col<3>() << nl;

        t3.col<0>({0, 2});
        Info<< "replaced col<0> = " << t3.col<0>() << nl;
        Info<< "tensor " << t3 << nl;

        t3.row<1>(Zero);
        Info<< "replaced row<1> = " << t3.row<1>() << nl;
        Info<< "tensor " << t3 << nl;
    }


    {
        vector2DField vfld1(8, Zero);

        forAll(vfld1, i)
        {
            vfld1[i] = (i+1) * ((i % 2) ? v1 : v2);
        }

        Info<< "vector: " << flatOutput(vfld1) << nl;

        scalarField xvals(8);
        scalarField yvals(8);
        unzip(vfld1, xvals, yvals);

        Info<< "unzip" << nl
            << " x => " << flatOutput(xvals) << nl
            << " y => " << flatOutput(yvals) << nl;

        reverse(xvals);
        zip(vfld1, xvals, yvals);

        Info<< "rezip (with reversed x)" << nl
            << "   => " << flatOutput(vfld1) << nl;
    }


    Info<< nl << "End\n" << nl;

    return 0;
}

// ************************************************************************* //
