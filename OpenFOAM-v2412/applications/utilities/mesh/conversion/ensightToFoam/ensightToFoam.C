/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    ensightToFoam

Group
    grpMeshConversionUtilities

Description
    Convert an Ensight Gold mesh into OpenFOAM format.

Usage
    \b ensightToFoam [OPTION] \<ensightGeometryFile\>

    Options:
      - \par -mergeTol \<factor\>
        Specify an alternative merging tolerance as a fraction of
        the bounding box of the points.

      - \par -scale \<factor\>
        Specify an optional geometry scaling factor.

      - \par -keepHandedness
        Do not automatically flip negative volume cells

See also
    Foam::meshReader and Foam::fileFormats::STARCDMeshReader

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "ensightMeshReader.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert Ensight mesh to OpenFOAM"
    );

    argList::noParallel();
    argList::addArgument(".geo file", "The file containing the geometry");
    argList::addOption
    (
        "mergeTol",
        "factor",
        "Merge tolerance as a fraction of bounding box - 0 to disable merging"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor - default is 1"
    );
    argList::addBoolOption
    (
        "keepHandedness",
        "Do not automatically flip inverted cells"
        " (default is to do a geometric test)"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    // More precision (for points data)
    IOstream::minPrecision(10);

    const fileName geomFile(args.get<fileName>(1));

    {
        fileFormats::ensightMeshReader reader
        (
            geomFile,
            runTime,
            args.getOrDefault<scalar>("mergeTol", 1e-10),
            args.getOrDefault<scalar>("scale", 1.0),
            args.found("keepHandedness")
        );

        autoPtr<polyMesh> mesh = reader.mesh(runTime);
        mesh().setInstance(runTime.constant());
        mesh().write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
