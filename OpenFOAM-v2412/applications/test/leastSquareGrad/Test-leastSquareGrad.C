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
    Test-leastSquareGrad

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "leastSquareGrad.H"
#include "labelVector.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    labelVector geomDim(1,1,-1);

    leastSquareGrad<scalar> lsGrad("polyDegree1", geomDim);

    List<vector> pos
    ({
        {0,0,0},
        {1,0,0},
        {-1,0,0},
        {0,1,0},
        {0,-1,0},
        // {1,1,0}
    });

    List<scalar> values
    ({
        0, 1, -1, 1, -1
        //, 2
    });


    Map<List<vector>> posMap;

    posMap.insert(0, pos);
    posMap.insert(1, pos);
    posMap.insert(2, pos);

    Map<List<scalar>> valMap;

    valMap.insert(0, values);
    valMap.insert(1, values);
    valMap.insert(2, values);

    Info<< lsGrad.grad(posMap, valMap) << nl;


    leastSquareGrad<vector> lsGradVec("polyDegree1", geomDim);

    List<vector> valuesVec
    ({
        {0,0,0},
        {1,0,0},
        {-1,0,0},
        {1,0,0},
        {-1,0,0}
    });

    Map<List<vector>> valMapVec;

    valMapVec.insert(0, valuesVec);
    valMapVec.insert(1, valuesVec);
    valMapVec.insert(2, valuesVec);

    Info<< lsGradVec.grad(posMap,valMapVec) << nl;

    return 0;
}


// ************************************************************************* //
