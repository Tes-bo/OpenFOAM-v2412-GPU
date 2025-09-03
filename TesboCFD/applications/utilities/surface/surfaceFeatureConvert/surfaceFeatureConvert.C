/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    surfaceFeatureConvert

Group
    grpSurfaceUtilities

Description
    Convert between edgeMesh formats.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "edgeMesh.H"

using namespace Foam;

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

    if (!edgeMesh::canReadType(readType, true))
    {
        good = false;
    }

    if (!edgeMesh::canWriteType(writeType, true))
    {
        good = false;
    }

    return good;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert between edgeMesh formats"
    );
    argList::noParallel();
    argList::addArgument("input", "The input edge file");
    argList::addArgument("output", "The output edge file");
    argList::addOption
    (
        "read-format",
        "type",
        "The input format (default: use file extension)"
    );
    argList::addOption
    (
        "write-format",
        "type",
        "The output format (default: use file extension)"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Input geometry scaling factor"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    const auto importName = args.get<fileName>(1);
    const auto exportName = args.get<fileName>(2);

    // Disable inplace editing
    if (importName == exportName)
    {
        FatalError
            << "Output file would overwrite input file."
            << exit(FatalError);
    }

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

    edgeMesh mesh(importName, readFileType);

    Info<< "\nRead edgeMesh " << importName << nl;
    mesh.writeStats(Info);
    Info<< nl
        << "\nwriting " << exportName;

    scalar scaleFactor(0);
    if (args.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< " with scaling " << scaleFactor << endl;
        mesh.scalePoints(scaleFactor);
    }
    else
    {
        Info<< " without scaling" << endl;
    }

    mesh.write(exportName, writeFileType);
    mesh.writeStats(Info);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
