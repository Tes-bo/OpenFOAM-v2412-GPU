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
    surfaceMeshImport

Group
    grpSurfaceUtilities

Description
    Import from various third-party surface formats into surfMesh
    with optional scaling or transformations (rotate/translate)
    on a coordinateSystem.

Usage
    \b surfaceMeshImport inputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface.

      - \par -name \<name\>
        Specify an alternative surface name when writing.

      - \par -read-format \<type\>
        Specify input file format

      - \par -read-scale \<scale\>
        Scale factor when reading files.

      - \par -write-scale \<scale\>
        Scale factor when writing files.

      - \par -dict \<dictionary\>
        Alternative dictionary for constant/coordinateSystems.

      - \par -from \<coordinateSystem\>
        Specify a coordinate system when reading files.

      - \par -to \<coordinateSystem\>
        Specify a coordinate system when writing files.

Note
    The filename extensions are used to determine the file format type.

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Import from various third-party surface formats into surfMesh"
    );

    argList::noParallel();
    argList::addArgument("surface", "The input surface file");

    argList::addBoolOption
    (
        "clean",
        "Perform some surface checking/cleanup on the input surface"
    );
    argList::addVerboseOption();
    argList::addOption
    (
        "name",
        "name",
        "The surface name when writing (default is 'default')"
    );
    argList::addOption
    (
        "read-format",
        "type",
        "Input format (default: use file extension)"
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

    #include "setRootCase.H"
    #include "createTime.H"

    // try for the latestTime, but create "constant" as needed
    instantList Times = runTime.times();
    if (Times.size())
    {
        label startTime = Times.size()-1;
        runTime.setTime(Times[startTime], startTime);
    }
    else
    {
        runTime.setTime(instant(0, runTime.constant()), 0);
    }


    const auto importName = args.get<fileName>(1);
    const auto exportName = args.getOrDefault<word>("name", "default");

    const int optVerbose  = args.verbose();

    const word readFileType
    (
        args.getOrDefault<word>("read-format", getExtension(importName))
    );

    // Check that reading is supported
    if (!meshedSurface::canRead(readFileType, true))
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
            << fromCsys->name() << endl;
        tmp<pointField> tpf = fromCsys->localPosition(surf.points());
        surf.movePoints(tpf());
    }

    if (toCsys)
    {
        Info<< "move points to coordinate system: "
            << toCsys->name() << endl;
        tmp<pointField> tpf = toCsys->globalPosition(surf.points());
        surf.movePoints(tpf());
    }

    if (args.readIfPresent("write-scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< "scale output " << scaleFactor << nl;
        surf.scalePoints(scaleFactor);
    }

    surfMesh smesh
    (
        IOobject
        (
            exportName,
            runTime.constant(),
            runTime
        ),
        std::move(surf)
    );


    Info<< "writing surfMesh:\n  " << smesh.objectPath() << endl;
    smesh.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
