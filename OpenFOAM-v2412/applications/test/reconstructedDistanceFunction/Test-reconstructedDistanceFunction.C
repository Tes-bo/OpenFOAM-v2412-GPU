/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |  www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 DLR
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "reconstructionSchemes.H"
#include "reconstructedDistanceFunction.H"
#include "Field.H"
#include "DynamicField.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field alpha1\n" << endl;

    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "create field phi\n" << endl;
    surfaceScalarField phi(fvc::interpolate(U) & mesh.Sf());

    dictionary dict = mesh.solverDict(alpha1.name());

    autoPtr<reconstructionSchemes> surf =
        reconstructionSchemes::New(alpha1,phi,U,dict);

    ++runTime;

    const volVectorField& centre = surf->centre();
    const volVectorField& normal = surf->normal();


    //pointField centres(0);
    //vectorField normals(0);
    DynamicField<point> centres(1000);
    DynamicField<vector> normals(1000);

    surf->reconstruct();

    zoneDistribute exchangeFields_(mesh);

    exchangeFields_.setUpCommforZone(surf->interfaceCell());

    Map<Field<vector>> mapCentres =
        exchangeFields_.getFields(surf->interfaceCell(),centre);

    Map<Field<vector>> mapNormal =
        exchangeFields_.getFields(surf->interfaceCell(),normal);

    forAll(surf->centre(),celli)
    {
        if (surf->interfaceCell()[celli])
        {
            centres.push_back(surf->centre()[celli]);
            normals.push_back(surf->normal()[celli]);
        }
    }

    reconstructedDistanceFunction distFunc(mesh);

    {
        runTime.cpuTimeIncrement();

        Info<< "Time " << runTime.cpuTimeIncrement()  << endl;

        #if 0
        distFunc.constructRDF
        (
            surf->interfaceCell(),
            surf->centre(),
            surf->normal(),
            2,
            exchangeFields_
        );
        #endif
    }

    runTime.write();

    return 0;
}


// ************************************************************************* //
