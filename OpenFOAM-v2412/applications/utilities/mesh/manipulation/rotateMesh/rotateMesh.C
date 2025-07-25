/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    rotateMesh

Group
    grpMeshManipulationUtilities

Description
    Rotates the mesh and fields from the direction n1 to direction n2.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "regionProperties.H"
#include "transformGeometricField.H"
#include "IOobjectList.H"

using namespace Foam;

template<class GeoField>
void ReadAndRotateFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const dimensionedTensor& rotT
)
{
    for (const IOobject& io : objects.csorted<GeoField>())
    {
        GeoField fld(io, mesh);
        Info<< "    Rotating " << fld.name() << endl;
        transform(fld, rotT, fld);
        fld.write();
    }
}


void rotateFields
(
    const fvMesh& mesh,
    const Time& runTime,
    const tensor& rotationT
)
{
    // Need dimensionedTensor for geometric fields
    const dimensionedTensor rotT(rotationT);

    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    ReadAndRotateFields<volVectorField>(mesh, objects, rotT);
    ReadAndRotateFields<volSphericalTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<volSymmTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<volTensorField>(mesh, objects, rotT);

    ReadAndRotateFields<surfaceVectorField>(mesh, objects, rotT);
    ReadAndRotateFields<surfaceSphericalTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<surfaceSymmTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<surfaceTensorField>(mesh, objects, rotT);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Rotate mesh points and vector/tensor fields\n"
        "Rotation from the <from> vector to the <to> vector"
    );

    timeSelector::addOptions();

    argList::addArgument("from", "The vector to rotate from");
    argList::addArgument("to",   "The vector to rotate to");

    #include "addAllRegionOptions.H"
    #include "setRootCase.H"

    const vector n1(args.get<vector>(1).normalise());
    const vector n2(args.get<vector>(2).normalise());

    const tensor rotT(rotationTensor(n1, n2));

    // ------------------------------------------------------------------------

    #include "createTime.H"

    // Handle -allRegions, -regions, -region
    #include "getAllRegionOptions.H"

    // ------------------------------------------------------------------------

    #include "createNamedMeshes.H"

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const fileName meshDir
        (
            polyMesh::meshDir(regionName)
        );

        if (regionNames.size() > 1)
        {
            Info<< "region=" << regionName << nl;
        }

        pointIOField points
        (
            IOobject
            (
                "points",
                runTime.findInstance(meshDir, "points"),
                meshDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        points = transform(rotT, points);

        // More precision (for points data)
        IOstream::minPrecision(10);

        Info<< "Writing points into directory "
            << runTime.relativePath(points.path()) << nl
            << endl;
        points.write();
    }

    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            if (regionNames.size() > 1)
            {
                Info<< "region=" << regionName << nl;
            }

            rotateFields(meshes[regioni], runTime, rotT);
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
