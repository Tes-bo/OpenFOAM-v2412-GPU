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
    surfaceSplitByPatch

Group
    grpSurfaceUtilities

Description
    Writes surface regions to separate files.

Usage
    \b surfaceSplitByPatch [OPTION]

    Options:
      - \par -patches NAME | LIST
        Specify single patch or multiple patches (name or regex) to extract
        For example,
        \verbatim
          -patches top
          -patches '( front \".*back\" )'
        \endverbatim

      - \par -exclude-patches NAME | LIST
        Exclude single or multiple patches (name or regex) from extracting.
        For example,
        \verbatim
          -exclude-patches '( inlet_1 inlet_2 "proc.*" )'
        \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "MeshedSurfaces.H"
#include "stringListOps.H"  // For stringListOps::findMatching()
#include "geometricSurfacePatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write surface mesh regions to separate files"
    );
    argList::noParallel();
    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to extract\n"
        "Eg, 'top' or '( front \".*back\" )'"
    );
    argList::addOption
    (
        "exclude-patches",
        "wordRes",
        "Exclude single or multiple patches (name or regex) from extracting.\n"
        "Eg, 'outlet' or '( inlet \".*Wall\" )'"
    );
    argList::addOptionCompat("exclude-patches", {"excludePatches", 2306});

    argList::addArgument("input", "The input surface file");

    argList args(argc, argv);

    const auto surfName = args.get<fileName>(1);

    const fileName surfBase(surfName.lessExt());

    const word extension(surfName.ext());


    Info<< nl
        << "Read surface from " << surfName << " ..." << nl << endl;

    meshedSurface surf(surfName);

    const surfZoneList& zones = surf.surfZones();

    Info<< "    " << surf.size() << " faces with "
        << zones.size() << " zones" << nl << nl;


    wordRes includePatches, excludePatches;

    if (args.readListIfPresent<wordRe>("patches", includePatches))
    {
        Info<< "Including patches " << flatOutput(includePatches)
            << nl << endl;
    }
    if (args.readListIfPresent<wordRe>("exclude-patches", excludePatches))
    {
        Info<< "Excluding patches " << flatOutput(excludePatches)
            << nl << endl;
    }

    // Identity if both include/exclude lists are empty
    const labelList zoneIndices
    (
        stringListOps::findMatching
        (
            zones,
            includePatches,
            excludePatches,
            nameOp<surfZone>()
        )
    );


    Info<< "Writing regions to "
        << zoneIndices.size() << " separate files ..." << nl << endl;


    // Faces to subset
    bitSet includeMap(surf.size());

    for (const label zonei : zoneIndices)
    {
        const surfZone& zn = zones[zonei];

        includeMap.reset();
        includeMap.set(zn.range());

        word patchName(zn.name());

        if (patchName.empty())
        {
            // In case people expect the same names as with triSurface
            patchName = geometricSurfacePatch::defaultName(zonei);
        }

        fileName outFile(surfBase + '_' + patchName + '.' + extension);

        Info<< "   Zone " << zonei << " (" << zn.size() << " faces) "
            << patchName
            << " to file " << outFile << nl;

        // Subset and write
        surf.subsetMesh(includeMap).write(outFile);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
