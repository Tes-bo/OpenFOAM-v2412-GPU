/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "NURBS3DVolumeCartesian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NURBS3DVolumeCartesian, 0);
    addToRunTimeSelectionTable
    (
        NURBS3DVolume,
        NURBS3DVolumeCartesian,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::NURBS3DVolumeCartesian::transformPointToCartesian
(
    const vector& localSystemCoordinates
) const
{
    return localSystemCoordinates;
}


Foam::tensor Foam::NURBS3DVolumeCartesian::transformationTensorDxDb
(
    label globalPointIndex
)
{
    return tensor::I;
}


void Foam::NURBS3DVolumeCartesian::updateLocalCoordinateSystem
(
    const vectorField& cartesianPoints
)
{
    localSystemCoordinates_ = cartesianPoints;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::NURBS3DVolumeCartesian::NURBS3DVolumeCartesian
(
    const dictionary& dict,
    const fvMesh& mesh,
    bool computeParamCoors
)
:
    NURBS3DVolume(dict, mesh, computeParamCoors)
{
    localSystemCoordinates_ = mesh_.points();
    writeCps("cpsBsplines" + mesh_.time().timeName());
    if (computeParamCoors)
    {
        getParametricCoordinates();
        determineActiveDesignVariablesAndPoints();
    }
}


// ************************************************************************* //
