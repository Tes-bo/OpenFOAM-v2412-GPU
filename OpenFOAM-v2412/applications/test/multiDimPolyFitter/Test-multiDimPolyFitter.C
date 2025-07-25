/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
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
    Test-multiDimPolyFitter

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiDimPolyFitter.H"
#include "labelVector.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    labelVector test(1,-1,-1);

    autoPtr<multiDimPolyFunctions> polyFunc
        = multiDimPolyFunctions::New("polyDegree1", test);

    vector vec(1,2,3);

    Info<< "polyFunc" << nl
        << "    nTerms " << polyFunc->nTerms() << nl
        << "    coeffs " << polyFunc->coeffs() << nl
        << "    termValues " << polyFunc->termValues(vec) << nl;

    multiDimPolyFitter<scalar> polyFitter("polyDegree1",test);

    List<vector> pos
    ({
        {0,0,0},
        {1,0,0},
        {-1,0,0}
    });

    List<scalar> values
    ({
        1, 2, 0
    });

    Info<< "pos " << pos << nl
        << "values " << values  << nl;

    scalarField fitData = polyFitter.fitData(pos, values);

    Info<< "fitData " << fitData  << nl;

    autoPtr<multiDimPolyFunctions> polyFuncDeg2
        = multiDimPolyFunctions::New("polyDegree2", labelVector(1,1,-1));

    //vector vec(1,2,3);

    Info<< "polyFunc" << nl
        << "    nTerms" << polyFuncDeg2->nTerms()  << nl
        // << "    coeffs " << polyFunc->coeffs()  << nl
        << "    termValues " << polyFuncDeg2->termValues(vector(1,1,0)) << nl;

    List<vector> pos2
    ({
        {0,0,0},
        {1,0,0},
        {-1,0,0},
        {0,1,0},
        {0,-1,0},
        {1,1,0}
    });

    List<scalar> values2
    ({
        0, 1, 1, 1, 1, 2
    });

    multiDimPolyFitter<scalar> polyFitter2
    {
        "polyDegree2", labelVector(1,1,-1)
    };

    scalarField fitData2 = polyFitter2.fitData(pos2, values2);

    Info<< "fitData " << fitData2  << nl;

    return 0;
}


// ************************************************************************* //
