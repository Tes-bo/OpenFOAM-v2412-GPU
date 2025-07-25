/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Test-surfaceWriter

Group
    grpSurfaceUtilities

Description
    Test surface writers.

Usage
    \b Test-surfaceWriter inputFile outputFile

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "surfaceWriter.H"
#include "MeshedSurfaces.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "convert between surface formats, "
        "but primarily for testing functionality\n"
        "Normally use surfaceMeshConvert instead."
    );

    argList::noParallel();
    argList::noFunctionObjects();

    argList::addOption
    (
        "type",
        "writerType"
    );

    argList::addArgument("inputFile");
    argList::addArgument("outputFile");

    #include "setRootCase.H"

    const auto importName = args.get<fileName>(1);
    const auto exportName = args.get<fileName>(2);

    if (importName == exportName)
    {
        FatalErrorInFunction
            << "Output file " << exportName << " would overwrite input file."
            << exit(FatalError);
    }


    if (!MeshedSurface<face>::canRead(importName, true))
    {
        return 1;
    }

    const word writerType = args.getOrDefault<word>("type", exportName.ext());

    auto surfWriter = surfaceWriter::New(writerType);

    {
        MeshedSurface<face> surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);

        Info<< "Open " << exportName
            << " for writing with " << surfWriter->type() << nl;

        surfWriter->open
        (
            surf.points(),
            surf.surfFaces(),
            exportName.lessExt(),
            false // serial
        );

        surfWriter->write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
