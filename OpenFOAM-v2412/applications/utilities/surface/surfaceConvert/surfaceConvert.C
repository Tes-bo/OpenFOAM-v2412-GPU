/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    surfaceConvert

Group
    grpSurfaceUtilities

Description
    Converts from one surface mesh format to another.

Usage
    \b surfaceConvert inputFile outputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface

      - \par -read-format \<type\>
        Specify input file format

      - \par -write-format \<type\>
        Specify output file format

      - \par -scale \<scale\>
        Specify a scaling factor for writing the files

      - \par -group
        Orders faces by region

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Time.H"

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

    if (!triSurface::canReadType(readType, true))
    {
        good = false;
    }

    if (!triSurface::canWriteType(writeType, true))
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
        "Convert between surface formats, using triSurface library components"
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("output", "The output surface file");

    argList::addBoolOption
    (
        "clean",
        "Perform some surface checking/cleanup on the input surface"
    );
    argList::addBoolOption
    (
        "group",
        "Reorder faces into groups; one per region"
    );
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
    argList::addOption
    (
        "precision",
        "int",
        "The output precision"
    );
    argList::addOptionCompat("precision", {"writePrecision", 1812});
    argList::addVerboseOption();

    argList args(argc, argv);

    {
        const unsigned prec = args.getOrDefault<unsigned>("precision", 0u);
        if (prec)
        {
            Info<< "Output write precision set to " << prec << endl;

            IOstream::defaultPrecision(prec);
            Sout.precision(prec);
        }
    }

    const auto importName = args.get<fileName>(1);
    const auto exportName = args.get<fileName>(2);

    const int optVerbose  = args.verbose();

    if (importName == exportName)
    {
        FatalError
            << "Output file would overwrite input file." << nl
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

    Info<< "Reading : " << importName << endl;
    triSurface surf(importName, readFileType, scaleFactor);

    if (args.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< "scale input " << scaleFactor << nl;
        surf.scalePoints(scaleFactor);
    }


    Info<< "Read surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    if (args.found("clean"))
    {
        Info<< "Cleaning up surface" << endl;
        surf.cleanup(optVerbose);

        Info<< "After cleaning up surface:" << endl;
        surf.writeStats(Info);
        Info<< endl;
    }

    const bool sortByRegion = args.found("group");
    if (sortByRegion)
    {
        Info<< "Reordering faces into groups; one per region." << endl;
    }
    else
    {
        Info<< "Maintaining face ordering" << endl;
    }

    Info<< "writing " << exportName << endl;

    surf.write(exportName, writeFileType, sortByRegion);

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
