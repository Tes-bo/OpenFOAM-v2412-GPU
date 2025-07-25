/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    surfaceMeshConvert

Group
    grpSurfaceUtilities

Description
    Convert between surface formats with optional scaling or
    transformations (rotate/translate) on a coordinateSystem.

Usage
    \b surfaceMeshConvert inputFile outputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface.

      - \par -read-format \<type\>
        The input file format (default: use file extension)

      - \par -write-format \<type\>
        The output file format (default: use file extension)

      - \par -read-scale \<scale\>
        Input geometry scaling factor.

      - \par -write-scale \<scale\>
        Output geometry scaling factor.

      - \par -dict \<dictionary\>
        Alternative dictionary for constant/coordinateSystems.

      - \par -from \<coordinateSystem\>
        Apply specified coordinate system after reading file.

      - \par -to \<coordinateSystem\>
        Apply specified coordinate system before writing file.

      - \par -tri
        Triangulate surface.

Note
    The filename extensions are used to determine the default file formats.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "MeshedSurfaces.H"
#include "coordinateSystems.H"
#include "cartesianCS.H"

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert between surface formats, using MeshSurface library components"
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("output", "The output surface file");

    argList::addBoolOption
    (
        "clean",
        "Perform some surface checking/cleanup on the input surface"
    );
    argList::addVerboseOption();
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
        "read-scale",
        "factor",
        "Input geometry scaling factor"
    );
    argList::addOption
    (
        "write-scale",
        "factor",
        "Output geometry scaling factor"
    );

    argList::addOptionCompat("read-scale", {"scaleIn", 1912});
    argList::addOptionCompat("write-scale", {"scaleOut", 1912});

    argList::addOption("dict", "file", "Alternative coordinateSystems");

    argList::addOption
    (
        "from",
        "system",
        "The source coordinate system, applied after '-read-scale'",
        true // advanced
    );
    argList::addOption
    (
        "to",
        "system",
        "The target coordinate system, applied before '-write-scale'",
        true // advanced
    );
    argList::addBoolOption
    (
        "tri",
        "Triangulate surface"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    const auto importName = args.get<fileName>(1);
    const auto exportName = args.get<fileName>(2);

    const int optVerbose  = args.verbose();

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


    scalar scaleFactor(0);

    // The coordinate transformations (must be cartesian)
    autoPtr<coordSystem::cartesian> fromCsys;
    autoPtr<coordSystem::cartesian> toCsys;

    if (args.found("from") || args.found("to"))
    {
        IOobject ioCsys = IOobject::selectIO
        (
            IOobject
            (
                coordinateSystems::typeName,
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            args.getOrDefault<fileName>("dict", "")
        );

        if (!ioCsys.typeHeaderOk<coordinateSystems>(false))
        {
            FatalError
                << "Cannot open coordinateSystems file\n    "
                << ioCsys.objectPath() << nl
                << exit(FatalError);
        }

        coordinateSystems globalCoords(ioCsys);

        if (args.found("from"))
        {
            const word csName(args["from"]);
            const auto* csPtr = globalCoords.cfind(csName);

            if (!csPtr)
            {
                FatalError
                    << "Cannot find -from " << csName << nl
                    << "available coordinateSystems: "
                    << flatOutput(globalCoords.names()) << nl
                    << exit(FatalError);
            }

            fromCsys = autoPtr<coordSystem::cartesian>::New(*csPtr);
        }

        if (args.found("to"))
        {
            const word csName(args["to"]);
            const auto* csPtr = globalCoords.cfind(csName);

            if (!csPtr)
            {
                FatalError
                    << "Cannot find -to " << csName << nl
                    << "available coordinateSystems: "
                    << flatOutput(globalCoords.names()) << nl
                    << exit(FatalError);
            }

            toCsys = autoPtr<coordSystem::cartesian>::New(*csPtr);
        }

        // Maybe fix this later
        if (fromCsys && toCsys)
        {
            FatalError
                << "Only allowed '-from' or '-to' option at the moment."
                << exit(FatalError);
        }
    }


    {
        meshedSurface surf(importName, readFileType);

        if (args.readIfPresent("read-scale", scaleFactor) && scaleFactor > 0)
        {
            Info<< "scale input " << scaleFactor << nl;
            surf.scalePoints(scaleFactor);
        }

        if (args.found("clean"))
        {
            surf.cleanup(optVerbose);
        }

        if (fromCsys)
        {
            Info<< "move points from coordinate system: "
                << fromCsys->name() << nl;
            tmp<pointField> tpf = fromCsys->localPosition(surf.points());
            surf.movePoints(tpf());
        }

        if (toCsys)
        {
            Info<< "move points to coordinate system: "
                << toCsys->name() << nl;
            tmp<pointField> tpf = toCsys->globalPosition(surf.points());
            surf.movePoints(tpf());
        }

        if (args.readIfPresent("write-scale", scaleFactor) && scaleFactor > 0)
        {
            Info<< "scale output " << scaleFactor << nl;
            surf.scalePoints(scaleFactor);
        }

        if (args.found("tri"))
        {
            Info<< "triangulate" << nl;
            surf.triangulate();
        }

        Info<< "writing " << exportName;
        surf.write(exportName, writeFileType);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
