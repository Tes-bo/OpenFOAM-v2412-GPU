/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022,2024 OpenCFD Ltd.
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
    surfaceTransformPoints

Group
    grpSurfaceUtilities

Description
    Transform (scale/rotate) a surface.
    Like transformPoints but for surfaces.

    The rollPitchYaw and yawPitchRoll options take three angles (degrees)
    that describe the intrinsic Euler rotation.

    rollPitchYaw
    - roll (rotation about X) followed by
    - pitch (rotation about Y) followed by
    - yaw (rotation about Z)

    yawPitchRoll
    - yaw (rotation about Z) followed by
    - pitch (rotation about Y) followed by
    - roll (rotation about X)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Fstream.H"
#include "boundBox.H"
#include "transformField.H"
#include "Pair.H"
#include "Tuple2.H"
#include "axisAngleRotation.H"
#include "EulerCoordinateRotation.H"
#include "MeshedSurfaces.H"
#include "cylindricalCS.H"

using namespace Foam;
using namespace Foam::coordinateRotations;

static word getExtension(const fileName& name)
{
    return
    (
        name.has_ext("gz")
      ? name.stem().ext()
      : name.ext()
    );
}


// Non-short-circuiting check to get all warnings
static bool hasReadWriteTypes(const word& readType, const word& writeType)
{
    volatile bool good = true;

    if (!meshedSurface::canReadType(readType, true))
    {
        good = false;
    }

    if (!meshedSurface::canWriteType(writeType, true))
    {
        good = false;
    }

    return good;
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
        "Transform (translate / rotate / scale) surface points.\n"
        "Like transformPoints but for surfaces.\n"
        "Note: roll=rotate about x, pitch=rotate about y, yaw=rotate about z"
    );
    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("output", "The output surface file");
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
    argList::addOption
    (
        "read-scale",
        "scalar | vector",
        "Uniform or non-uniform input scaling"
    );
    argList::addOption
    (
        "write-scale",
        "scalar | vector",
        "Uniform or non-uniform output scaling"
    );
    argList::addOption
    (
        "read-format",
        "type",
        "Input format (default: use file extension)"
    );
    argList::addOption
    (
        "write-format",
        "type",
        "Output format (default: use file extension)"
    );
    argList::addOption
    (
        "cylToCart",
        "(originVec axisVec directionVec)",
        "Tranform cylindrical coordinates to cartesian coordinates"
    );

    // Backward compatibility and with transformPoints
    argList::addOptionCompat("write-scale", {"scale", -2006});

    argList args(argc, argv);

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
            "read-scale",
            "write-scale",
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

    const auto importName = args.get<fileName>(1);
    const auto exportName = args.get<fileName>(2);

    const word readFileType
    (
        args.getOrDefault<word>("read-format", getExtension(importName))
    );

    const word writeFileType
    (
        args.getOrDefault<word>("write-format", getExtension(exportName))
    );


    // Check that reading/writing is supported
    if (!hasReadWriteTypes(readFileType, writeFileType))
    {
        FatalError
            << "Unsupported file format(s)" << nl
            << exit(FatalError);
    }


    Info<< "Reading surf from " << importName << " ..." << nl
        << "Writing surf to " << exportName << " ..." << endl;


    meshedSurface surf1(importName, readFileType);

    pointField points(surf1.points());


    // Begin operations

    // Input scaling
    applyScaling(points, getScalingOpt("read-scale", args));

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
    }

    if (useRotationCentre)
    {
        Info<< "Unset centre of rotation from " << rotationCentre << endl;
        points += rotationCentre;
    }

    // Output scaling
    applyScaling(points, getScalingOpt("write-scale", args));

    // Conversion to cylindrical coords
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

    surf1.movePoints(points);
    surf1.write(exportName, writeFileType);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
