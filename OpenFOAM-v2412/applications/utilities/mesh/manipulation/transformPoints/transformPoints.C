/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
    Copyright (C) 2024 Haakan Nilsson
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
    transformPoints

Group
    grpMeshManipulationUtilities

Description
    Transforms the mesh points in the polyMesh directory according to the
    translate, rotate and scale options.

Usage
    Options are:

    -time value
        Specify the time to search from and apply the transformation
        (default is latest)

    -recentre
        Recentre using the bounding box centre before other operations

    -translate vector
        Translate the points by the given vector before rotations

    -rotate (vector vector)
        Rotate the points from the first vector to the second

    -rotate-angle (vector angle)
        Rotate angle degrees about vector axis.

    -rotate-x angle
        Rotate (degrees) about x-axis.

    -rotate-y angle
        Rotate (degrees) about y-axis.

    -rotate-z angle
        Rotate (degrees) about z-axis.

     or -yawPitchRoll : (yaw pitch roll) degrees
     or -rollPitchYaw : (roll pitch yaw) degrees

    -scale scalar|vector
        Scale the points by the given scalar or vector on output.

    The any or all of the three options may be specified and are processed
    in the above order.

    -cylToCart (originVector axisVector directionVector)
        Tranform cylindrical coordinates to cartesian coordinates

    With -rotateFields (in combination with -rotate/yawPitchRoll/rollPitchYaw)
    it will also read & transform vector & tensor fields.

Note
    roll (rotation about x)
    pitch (rotation about y)
    yaw (rotation about z)

    - with -rotate and two exactly opposing vectors it will actually mirror
    the geometry. Use any of the other rotation options instead or use
    two steps, each 90 degrees.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"
#include "regionProperties.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "axisAngleRotation.H"
#include "EulerCoordinateRotation.H"
#include "cylindricalCS.H"

using namespace Foam;
using namespace Foam::coordinateRotations;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void readAndRotateFields
(
    PtrList<GeoField>& flds,
    const fvMesh& mesh,
    const dimensionedTensor& rotT,
    const IOobjectList& objects
)
{
    ReadFields(mesh, objects, flds);
    for (GeoField& fld : flds)
    {
        Info<< "Transforming " << fld.name() << endl;
        transform(fld, rotT, fld);
    }
}


void rotateFields
(
    const word& regionName,
    const Time& runTime,
    const tensor& rotT
)
{
    // Need dimensionedTensor for geometric fields
    const dimensionedTensor T(rotT);

    #include "createRegionMesh.H"

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    readAndRotateFields(vsFlds, mesh, T, objects);

    PtrList<volVectorField> vvFlds;
    readAndRotateFields(vvFlds, mesh, T, objects);

    PtrList<volSphericalTensorField> vstFlds;
    readAndRotateFields(vstFlds, mesh, T, objects);

    PtrList<volSymmTensorField> vsymtFlds;
    readAndRotateFields(vsymtFlds, mesh, T, objects);

    PtrList<volTensorField> vtFlds;
    readAndRotateFields(vtFlds, mesh, T, objects);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    readAndRotateFields(ssFlds, mesh, T, objects);

    PtrList<surfaceVectorField> svFlds;
    readAndRotateFields(svFlds, mesh, T, objects);

    PtrList<surfaceSphericalTensorField> sstFlds;
    readAndRotateFields(sstFlds, mesh, T, objects);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    readAndRotateFields(ssymtFlds, mesh, T, objects);

    PtrList<surfaceTensorField> stFlds;
    readAndRotateFields(stFlds, mesh, T, objects);

    mesh.write();
}


// Retrieve scaling option
// - size 0 : no scaling
// - size 1 : uniform scaling
// - size 3 : non-uniform scaling
List<scalar> getScalingOpt(const word& optName, const argList& args)
{
    // readListIfPresent handles single or multiple values
    // - accept 1 or 3 values

    List<scalar> scaling;
    args.readListIfPresent(optName, scaling);

    if (scaling.size() == 1)
    {
        // Uniform scaling
    }
    else if (scaling.size() == 3)
    {
        // Non-uniform, but may actually be uniform
        if
        (
            equal(scaling[0], scaling[1])
         && equal(scaling[0], scaling[2])
        )
        {
            scaling.resize(1);
        }
    }
    else if (!scaling.empty())
    {
        FatalError
            << "Incorrect number of components, must be 1 or 3." << nl
            << "    -" << optName << ' ' << args[optName].c_str() << endl
            << exit(FatalError);
    }

    if (scaling.size() == 1 && equal(scaling[0], 1))
    {
        // Scale factor 1 == no scaling
        scaling.clear();
    }

    // Zero and negative scaling are permitted

    return scaling;
}


void applyScaling(pointField& points, const List<scalar>& scaling)
{
    if (scaling.size() == 1)
    {
        Info<< "Scaling points uniformly by " << scaling[0] << nl;
        points *= scaling[0];
    }
    else if (scaling.size() == 3)
    {
        const vector factor(scaling[0], scaling[1], scaling[2]);
        Info<< "Scaling points by " << factor << nl;
        cmptMultiply(points, points, factor);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transform (translate / rotate / scale) mesh points.\n"
        "Note: roll=rotate about x, pitch=rotate about y, yaw=rotate about z"
    );
    argList::addOption
    (
        "time",
        "time",
        "Specify the time to search from and apply the transformation"
        " (default is latest)"
    );
    argList::addBoolOption
    (
        "recentre",
        "Recentre the bounding box before other operations"
    );
    argList::addOption
    (
        "translate",
        "vector",
        "Translate by specified <vector> before rotations"
    );
    argList::addBoolOption
    (
        "auto-centre",
        "Use bounding box centre as centre for rotations"
    );
    argList::addOption
    (
        "centre",
        "point",
        "Use specified <point> as centre for rotations"
    );
    argList::addOptionCompat("auto-centre", {"auto-origin", 2206});
    argList::addOptionCompat("centre", {"origin", 2206});

    argList::addOption
    (
        "rotate",
        "(vectorA vectorB)",
        "Rotate from <vectorA> to <vectorB> - eg, '((1 0 0) (0 0 1))'"
    );
    argList::addOption
    (
        "rotate-angle",
        "(vector angle)",
        "Rotate <angle> degrees about <vector> - eg, '((1 0 0) 45)'"
    );
    argList::addOption
    (
        "rotate-x", "deg",
        "Rotate (degrees) about x-axis"
    );
    argList::addOption
    (
        "rotate-y", "deg",
        "Rotate (degrees) about y-axis"
    );
    argList::addOption
    (
        "rotate-z", "deg",
        "Rotate (degrees) about z-axis"
    );
    argList::addOption
    (
        "rollPitchYaw",
        "vector",
        "Rotate by '(roll pitch yaw)' degrees"
    );
    argList::addOption
    (
        "yawPitchRoll",
        "vector",
        "Rotate by '(yaw pitch roll)' degrees"
    );
    argList::addBoolOption
    (
        "rotateFields",
        "Read and transform vector and tensor fields too"
    );
    argList::addOption
    (
        "scale",
        "scalar | vector",
        "Scale by the specified amount - Eg, for uniform [mm] to [m] scaling "
        "use either '(0.001 0.001 0.001)' or simply '0.001'"
    );
    argList::addOption
    (
        "cylToCart",
        "(originVec axisVec directionVec)",
        "Tranform cylindrical coordinates to cartesian coordinates"
    );


    // Compatibility with surfaceTransformPoints
    argList::addOptionCompat("scale", {"write-scale", 0});

    #include "addAllRegionOptions.H"
    #include "setRootCase.H"

    const bool doRotateFields = args.found("rotateFields");

    // Verify that an operation has been specified
    {
        const List<word> operationNames
        ({
            "recentre",
            "translate",
            "rotate",
            "rotate-angle",
            "rotate-x",
            "rotate-y",
            "rotate-z",
            "rollPitchYaw",
            "yawPitchRoll",
            "scale",
            "cylToCart"
        });

        if (!args.count(operationNames))
        {
            FatalError
                << "No operation supplied, "
                << "use at least one of the following:" << nl
                << "   ";

            for (const auto& opName : operationNames)
            {
                FatalError
                    << " -" << opName;
            }

            FatalError
                << nl << exit(FatalError);
        }
    }

    // ------------------------------------------------------------------------

    #include "createTime.H"

    if (args.found("time"))
    {
        if (args["time"] == "constant")
        {
            runTime.setTime(instant(0, "constant"), 0);
        }
        else
        {
            const scalar timeValue = args.get<scalar>("time");
            runTime.setTime(instant(timeValue), 0);
        }
    }

    // Handle -allRegions, -regions, -region
    #include "getAllRegionOptions.H"

    // ------------------------------------------------------------------------

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


        // Begin operations

        vector v;
        if (args.found("recentre"))
        {
            v = boundBox(points).centre();
            Info<< "Adjust centre " << v << " -> (0 0 0)" << endl;
            points -= v;
        }

        if (args.readIfPresent("translate", v))
        {
            Info<< "Translating points by " << v << endl;
            points += v;
        }

        vector rotationCentre;
        bool useRotationCentre = args.readIfPresent("centre", rotationCentre);
        if (args.found("auto-centre") && !useRotationCentre)
        {
            useRotationCentre = true;
            rotationCentre = boundBox(points).centre();
        }

        if (useRotationCentre)
        {
            Info<< "Set centre of rotation to " << rotationCentre << endl;
            points -= rotationCentre;
        }


        // Get a rotation specification

        tensor rot(Zero);
        bool useRotation(false);

        if (args.found("rotate"))
        {
            Pair<vector> n1n2
            (
                args.lookup("rotate")()
            );
            n1n2[0].normalise();
            n1n2[1].normalise();

            rot = rotationTensor(n1n2[0], n1n2[1]);
            useRotation = true;
        }
        else if (args.found("rotate-angle"))
        {
            const Tuple2<vector, scalar> rotAxisAngle
            (
                args.lookup("rotate-angle")()
            );

            const vector& axis = rotAxisAngle.first();
            const scalar angle = rotAxisAngle.second();

            Info<< "Rotating points " << nl
                << "    about " << axis << nl
                << "    angle " << angle << nl;

            rot = axisAngle::rotation(axis, angle, true);
            useRotation = true;
        }
        else if (args.found("rotate-x"))
        {
            const scalar angle = args.get<scalar>("rotate-x");

            Info<< "Rotating points about x-axis: " << angle << nl;

            rot = axisAngle::rotation(vector::X, angle, true);
            useRotation = true;
        }
        else if (args.found("rotate-y"))
        {
            const scalar angle = args.get<scalar>("rotate-y");

            Info<< "Rotating points about y-axis: " << angle << nl;

            rot = axisAngle::rotation(vector::Y, angle, true);
            useRotation = true;
        }
        else if (args.found("rotate-z"))
        {
            const scalar angle = args.get<scalar>("rotate-z");

            Info<< "Rotating points about z-axis: " << angle << nl;

            rot = axisAngle::rotation(vector::Z, angle, true);
            useRotation = true;
        }
        else if (args.readIfPresent("rollPitchYaw", v))
        {
            Info<< "Rotating points by" << nl
                << "    roll  " << v.x() << nl
                << "    pitch " << v.y() << nl
                << "    yaw   " << v.z() << nl;

            rot = euler::rotation(euler::eulerOrder::ROLL_PITCH_YAW, v, true);
            useRotation = true;
        }
        else if (args.readIfPresent("yawPitchRoll", v))
        {
            Info<< "Rotating points by" << nl
                << "    yaw   " << v.x() << nl
                << "    pitch " << v.y() << nl
                << "    roll  " << v.z() << nl;

            rot = euler::rotation(euler::eulerOrder::YAW_PITCH_ROLL, v, true);
            useRotation = true;
        }

        if (useRotation)
        {
            Info<< "Rotating points by " << rot << endl;
            transform(points, rot, points);

            if (doRotateFields)
            {
                rotateFields(regionName, runTime, rot);
            }
        }

        if (useRotationCentre)
        {
            Info<< "Unset centre of rotation from " << rotationCentre << endl;
            points += rotationCentre;
        }

        // Output scaling
        applyScaling(points, getScalingOpt("scale", args));

        if (args.found("cylToCart"))
        {
            vectorField n1n2(args.lookup("cylToCart")());
            n1n2[1].normalise();
            n1n2[2].normalise();

            cylindricalCS ccs
                (
                    "ccs",
                    n1n2[0],
                    n1n2[1],
                    n1n2[2]
                );

            points = ccs.globalPosition(points);
        }

        // More precision (for points data)
        IOstream::minPrecision(10);

        Info<< "Writing points into directory "
            << runTime.relativePath(points.path()) << nl
            << endl;
        points.write();
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
