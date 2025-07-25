/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Test-quaternion

Description
    Test application for quaternions.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "quaternion.H"
#include "septernion.H"
#include "unitConversion.H"
#include "Tuple2.H"
#include "IOstreams.H"
#include "transform.H"
#include "axisAngleRotation.H"
#include "EulerCoordinateRotation.H"

using namespace Foam;
using namespace Foam::coordinateRotations;


void printRotation(const tensor& rot)
{
    Info<< "[\n"
        << "    " << rot.xx() << ' ' << rot.xy() << ' ' << rot.xz() << nl
        << "    " << rot.yx() << ' ' << rot.yy() << ' ' << rot.yz() << nl
        << "    " << rot.zx() << ' ' << rot.zy() << ' ' << rot.zz() << nl
        << "]\n";
}


void printRotation(const quaternion& quat)
{
    tensor rot(quat.R());

    Info<< "quaternion " << quat << nl
        << "rotation" << nl;

    printRotation(rot);
    Info<< "transpose" << nl;
    printRotation(rot.T());
}


bool equalTensors(const tensor& rot1, const tensor& rot2)
{
    for (direction cmpt=0; cmpt < tensor::nComponents; ++cmpt)
    {
        // Cannot be really picky, but SMALL is reasonable
        if (mag(rot1[cmpt] - rot2[cmpt]) > SMALL)
        {
            return false;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "rollPitchYaw",
        "vector",
        "Rotate by '(roll pitch yaw)' in degrees"
    );
    argList::addOption
    (
        "yawPitchRoll",
        "vector",
        "Rotate by '(yaw pitch roll)' in degrees"
    );
    argList::addOption
    (
        "euler",
        "vector",
        "Rotate by '(phi theta psi)' in degrees"
    );
    argList::addOption
    (
        "xyz",
        "vector",
        "Rotate about x-y-z axes in degrees"
    );
    argList::addOption
    (
        "zyx",
        "vector",
        "Rotate about z-y-x axes in degrees"
    );
    argList::addOption
    (
        "rotate-angle",
        "(vector angle)",
        "Rotate about the <vector> by <angle> degrees - eg, '((1 0 0) 45)'"
    );
    argList::addVerboseOption
    (
        "Report euler angles"
    );

    argList args(argc, argv);


    vector rotVector;

    if (args.readIfPresent("rollPitchYaw", rotVector))
    {
        Info<< nl
            << "Rotate by" << nl
            << "    roll  " << rotVector.x() << nl
            << "    pitch " << rotVector.y() << nl
            << "    yaw   " << rotVector.z() << nl;

        rotVector *= degToRad();

        const quaternion quat(quaternion::eulerOrder::XYZ, rotVector);

        printRotation(quat);

        // Euler
        const tensor rot
        (
            euler::rotation(euler::eulerOrder::XYZ, rotVector, false)
        );
        printRotation(rot);
    }

    if (args.readIfPresent("yawPitchRoll", rotVector))
    {
        Info<< nl
            << "Rotate by" << nl
            << "    yaw   " << rotVector.x() << nl
            << "    pitch " << rotVector.y() << nl
            << "    roll  " << rotVector.z() << nl;

        rotVector *= degToRad();

        const quaternion quat(quaternion::eulerOrder::ZYX, rotVector);

        printRotation(quat);

        // Euler
        const tensor rot
        (
            euler::rotation(euler::eulerOrder::ZYX, rotVector, false)
        );
        printRotation(rot);
    }
    if (args.readIfPresent("euler", rotVector))
    {
        Info<< nl
            << "Rotate by" << nl
            << "    phi   " << rotVector.x() << nl
            << "    theta " << rotVector.y() << nl
            << "    psi   " << rotVector.z() << nl;

        printRotation(euler(rotVector, true).R());

        rotVector *= degToRad();

        const quaternion quat(quaternion::eulerOrder::ZXZ, rotVector);

        printRotation(quat);
    }
    if (args.readIfPresent("xyz", rotVector))
    {
        Info<< nl
            << "Rotate about" << nl
            << "    x   " << rotVector.x() << nl
            << "    y   " << rotVector.y() << nl
            << "    z   " << rotVector.z() << nl;

        Vector<tensor> Rs
        (
            axisAngle(vector(1,0,0), rotVector.x(), true).R(),
            axisAngle(vector(0,1,0), rotVector.y(), true).R(),
            axisAngle(vector(0,0,1), rotVector.z(), true).R()
        );

        printRotation(Rs.x() & Rs.y() & Rs.z());
    }
    if (args.readIfPresent("zyx", rotVector))
    {
        Info<< nl
            << "Rotate about" << nl
            << "    z   " << rotVector.x() << nl
            << "    y   " << rotVector.y() << nl
            << "    x   " << rotVector.z() << nl;

        Vector<tensor> Rs
        (
            axisAngle(vector(1,0,0), rotVector.z(), true).R(),
            axisAngle(vector(0,1,0), rotVector.y(), true).R(),
            axisAngle(vector(0,0,1), rotVector.x(), true).R()
        );

        printRotation(Rs.z() & Rs.y() & Rs.x());
    }

    if (args.found("rotate-angle"))
    {
        const Tuple2<vector, scalar> axisAngle
        (
            args.lookup("rotate-angle")()
        );

        Info<< nl
            << "Rotate" << nl
            << "    about " << axisAngle.first() << nl
            << "    angle " << axisAngle.second() << nl;

        const quaternion quat
        (
            axisAngle.first(),
            degToRad(axisAngle.second())
        );

        printRotation(quat);

        Info<< "transform Ra = "
            << Ra
               (
                   axisAngle.first() / mag(axisAngle.first()),
                   degToRad(axisAngle.second())
               ) << nl;
        Info<< "-ve Ra = "
            << Ra
               (
                   axisAngle.first() / mag(axisAngle.first()),
                   degToRad(-axisAngle.second())
               ) << nl;
    }


    Info<< nl << nl;

    quaternion q(vector(1, 2, 3), 0.7853981);
    Info<< "q " << q << nl;

    vector v(0.1, 0.4, 2.1);
    Info<< "v " << v << nl;

    Info<< "inv(q)*q " << inv(q)*q << nl;

    Info<< "q*quaternion(0, v)*conjugate(q) "
        << q*quaternion(0, v)*conjugate(q) << nl;

    Info<< "q.R() " << q.R() << nl;
    Info<< "q.transform(v) " << q.transform(v) << nl;
    Info<< "q.R() & v " << (q.R() & v) << nl;
    Info<< "quaternion(q.R()).transform(v) "
        << (quaternion(q.R()).transform(v)) << nl;

    Info<< "q.invTransform(v) " << q.invTransform(v) << nl;

    septernion tr(vector(0, 0.1, 0), q);
    Info<< "tr " << tr << nl;

    Info<< "inv(tr)*tr " << inv(tr)*tr << nl;

    Info<< "tr.transform(v) " << tr.transformPoint(v) << nl;

    vector origin(1, 2, 4);

    Info<< "(septernion(-origin)*q*septernion(origin))"
        << ".transform(v) "
        << (septernion(-origin)*q*septernion(origin)).transformPoint(v)
        <<  " "
        << septernion(-origin)
          .transformPoint(q.transform(septernion(origin).transformPoint(v)))
        << nl;

    Info<< "Test conversion from and to Euler-angles" << nl;
    vector angles(0.1, 0.2, 0.3);

    for (int rs : quaternion::eulerOrderNames.values())
    {
        const quaternion::eulerOrder order = quaternion::eulerOrder(rs);
        const word& orderName = quaternion::eulerOrderNames[order];

        if
        (
            mag(angles - quaternion(order, angles).eulerAngles(order))
          > SMALL
        )
        {
            FatalErrorInFunction
                << "Inconsistent conversion for euler rotation order "
                << orderName << nl << exit(FatalError);
        }


        tensor rotQ(quaternion(order, angles).R());
        tensor rotE(euler::rotation(order, angles, false));

        if (args.verbose())
        {
            Info<< "euler " << orderName << angles << nl;
            printRotation(rotE);
        }

        if (!equalTensors(rotQ, rotE))
        {
            WarningInFunction
                << "Inconsistent quaternion/euler rotation matrices for "
                << orderName << nl;

            printRotation(rotQ);
            printRotation(rotE);

            FatalError
                << nl << exit(FatalError);
        }
    }

    List<septernion> ss(3);
    List<scalar> w(3);

    ss[0] = septernion(vector(0, 0.1, 0), quaternion(0.7, vector(1, 2, 3)));
    w[0] = 0.1;
    ss[1] = septernion(vector(0, 0.2, 0), quaternion(-0.6, vector(-2, -1, -3)));
    w[1] = 0.5;
    ss[2] = septernion(vector(0, 0.3, 0), quaternion(0.3, vector(3, 2, 1)));
    w[2] = 0.4;

    Info<< "average(ss, w) " << average(ss, w) << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
