/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    mapDistributePolyMesh

Description
    Test for procAddressing

\*---------------------------------------------------------------------------*/

#include "IOmapDistributePolyMesh.H"
#include "argList.H"
#include "Time.H"
#include "surfaceFields.H"
#include "flipOp.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    Info<< "Reading distribute map\n" << endl;

    const word instance("0.005");
    const scalar instanceValue(0.005);


    IOmapDistributePolyMesh map
    (
        IOobject
        (
            "procAddressing",
            instance,
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    {
        // Load the instance mesh
        runTime.setTime(instanceValue, 0);
        polyMesh distributedMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                instance,
                runTime,
                IOobject::MUST_READ
            )
        );

        // Faces, no flip
        {
            const mapDistribute& faceMap = map.faceMap();
            pointField fc(mesh.faceCentres());
            faceMap.distribute(fc, identityOp());
            Pout<< "Construct size:" << faceMap.constructSize() << endl;
            forAll(distributedMesh.faceCentres(), facei)
            {
                Pout<< "face:" << facei
                    << "\tmappedFc:" << fc[facei]
                    << "\tactual:" << distributedMesh.faceCentres()[facei]
                    << endl;
            }
        }
        // Faces, flipped field
        {
            const mapDistribute& faceMap = map.faceMap();
            scalarField flux(mesh.faceAreas() & vector(1, 1, 1));
            faceMap.distribute(flux, flipOp());
            Pout<< "Construct size:" << faceMap.constructSize() << endl;
            const scalarField newFlux
            (
                distributedMesh.faceAreas()
              & vector(1, 1, 1)
            );
            forAll(newFlux, facei)
            {
                Pout<< "face:" << facei
                    << "\tmappedFlux:" << flux[facei]
                    << "\tactual:" << newFlux[facei]
                    << endl;
            }
        }

        {
            const mapDistribute& cellMap = map.cellMap();
            pointField cc(mesh.cellCentres());
            cellMap.distribute(cc, identityOp());
            Pout<< "Construct size:" << cellMap.constructSize() << endl;
            forAll(distributedMesh.cellCentres(), celli)
            {
                Pout<< "cell:" << celli
                    << "\tmappedCc:" << cc[celli]
                    << "\tactual:" << distributedMesh.cellCentres()[celli]
                    << endl;
            }
        }
        {
            const mapDistribute& pointMap = map.pointMap();
            pointField pc(mesh.points());
            pointMap.distribute(pc, identityOp());
            Pout<< "Construct size:" << pointMap.constructSize() << endl;
            forAll(distributedMesh.points(), pointi)
            {
                Pout<< "point:" << pointi
                    << "\tmappedPoint:" << pc[pointi]
                    << "\tactual:" << distributedMesh.points()[pointi]
                    << endl;
            }
        }
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
