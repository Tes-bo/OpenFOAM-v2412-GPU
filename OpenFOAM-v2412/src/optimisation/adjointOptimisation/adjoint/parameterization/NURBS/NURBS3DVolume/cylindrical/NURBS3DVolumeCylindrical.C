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

#include "NURBS3DVolumeCylindrical.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointPatchFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(NURBS3DVolumeCylindrical, 0);
    addToRunTimeSelectionTable
    (
        NURBS3DVolume,
        NURBS3DVolumeCylindrical,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::NURBS3DVolumeCylindrical::transformPointToCartesian
(
    const vector& localSystemCoordinates
) const
{
    vector cartesianCoors
    (
        localSystemCoordinates.x()*cos(localSystemCoordinates.y()),
        localSystemCoordinates.x()*sin(localSystemCoordinates.y()),
        localSystemCoordinates.z()
    );
    cartesianCoors += origin_;

    return cartesianCoors;
}


Foam::tensor Foam::NURBS3DVolumeCylindrical::transformationTensorDxDb
(
    label globalPointIndex
)
{
    const vector& localCoors = localSystemCoordinates_[globalPointIndex];
    tensor transformTensor
    (
        cos(localCoors.y()), -localCoors.x()*sin(localCoors.y()), 0,
        sin(localCoors.y()), localCoors.x()*cos(localCoors.y()), 0,
        0, 0, 1
    );

    return transformTensor;
}


void Foam::NURBS3DVolumeCylindrical::updateLocalCoordinateSystem
(
    const vectorField& cartesianPoints
)
{
    forAll(cartesianPoints, pI)
    {
        const vector point(cartesianPoints[pI] - origin_);
        vector cylindricalCoors(Zero);

        const scalar R(Foam::sqrt(sqr(point.x()) + sqr(point.y())));
        const scalar theta(atan2(point.y(), point.x()));
        cylindricalCoors.x() = R;
        cylindricalCoors.y() = theta;
        cylindricalCoors.z() = cartesianPoints[pI].z();
        localSystemCoordinates_[pI] = cylindricalCoors;
    }

    pointVectorField cylindricalCoors
    (
        IOobject
        (
           "cylindricalCoors" + name_,
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,
           IOobject::NO_WRITE
        ),
        pointMesh::New(mesh_),
        vector::zero
    );
    cylindricalCoors.primitiveFieldRef() = localSystemCoordinates_;
    cylindricalCoors.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NURBS3DVolumeCylindrical::NURBS3DVolumeCylindrical
(
    const dictionary& dict,
    const fvMesh& mesh,
    bool computeParamCoors
)
:
    NURBS3DVolume(dict, mesh, computeParamCoors),
    origin_(dict.get<vector>("origin"))
{
    updateLocalCoordinateSystem(mesh.points());
    writeCps("cpsBsplines" + mesh_.time().timeName());
    if (computeParamCoors)
    {
        getParametricCoordinates();
        determineActiveDesignVariablesAndPoints();
    }
}


// ************************************************************************* //
