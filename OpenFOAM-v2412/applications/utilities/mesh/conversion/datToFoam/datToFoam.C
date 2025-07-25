/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    datToFoam

Group
    grpMeshConversionUtilities

Description
    Reads in a datToFoam mesh file and outputs a points file.
    Used in conjunction with blockMesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Fstream.H"
#include "polyMesh.H"
#include "unitConversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reads in a datToFoam mesh file and outputs a points file.\n"
        "Used in conjunction with blockMesh."
    );

    argList::noParallel();
    argList::addArgument("dat file");

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    #include "createTime.H"

    std::ifstream plot3dFile(args.get<fileName>(1));

    string line;
    std::getline(plot3dFile, line);
    std::getline(plot3dFile, line);

    IStringStream Istring(line);
    word block;
    string zoneName;
    token punctuation;
    label iPoints, jPoints;

    Istring >> block;
    Istring >> block;
    Istring >> zoneName;
    Istring >> punctuation;
    Istring >> block;
    Istring >> iPoints;
    Istring >> block;
    Istring >> jPoints;

    Info<< "Number of vertices in i direction = " << iPoints << nl
        << "Number of vertices in j direction = " << jPoints << endl;

    // We ignore the first layer of points in i and j the biconic meshes
    const label nPointsij = (iPoints - 1)*(jPoints - 1);

    pointField points(nPointsij, Zero);

    for (direction comp = 0; comp < 2; ++comp)
    {
        label p = 0;

        for (label j = 0; j < jPoints; ++j)
        {
            for (label i = 0; i < iPoints; ++i)
            {
                scalar coord;
                plot3dFile >> coord;

                // if statement ignores the first layer in i and j
                if (i>0 && j>0)
                {
                    points[p++][comp] = coord;
                }
            }
        }
    }

    // correct error in biconic meshes
    for (point& p : points)
    {
        if (p.y() < 1e-07)
        {
            p.y() = 0;
        }
    }

    pointField pointsWedge(nPointsij*2, Zero);

    const scalar a = degToRad(0.1);
    const tensor rotateZ
    (
        1.0,  0.0,      0.0,
        0.0,  1.0,      0.0,
        0.0, -::sin(a), ::cos(a)
    );

    forAll(points, i)
    {
        pointsWedge[i] = (rotateZ & points[i]);
        pointsWedge[i + nPointsij] =
            cmptMultiply
            (
                vector(1.0, 1.0, -1.0),
                pointsWedge[i]
            );
    }

    // Write (overwrite) points
    {
        IOobject iopoints
        (
            "points",
            runTime.constant(),     // instance
            polyMesh::meshSubDir,   // local
            runTime
        );

        if (!exists(iopoints.path()))
        {
            mkDir(iopoints.path());
        }

        OFstream os(iopoints.objectPath(), runTime.writeStreamOption());

        Info<< "Writing points to: " << nl
            << "    " << os.name() << endl;

        iopoints.writeHeader(os, vectorIOField::typeName);
        os << pointsWedge;
        IOobject::writeEndDivider(os);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
