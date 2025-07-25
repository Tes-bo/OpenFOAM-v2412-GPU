/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    ccmToFoam

Group
    grpMeshConversionUtilities

Description
    Reads CCM files as written by PROSTAR/STARCCM and writes an
    OPENFOAM polyMesh.

Usage
    \b ccmToFoam [OPTION] ccmMesh

    Options:
      - \par -ascii
        Write in ASCII format instead of binary

      - \par -export
        re-export mesh in CCM format for post-processing

      - \par -list
        List some information about the geometry

      - \par -name \<name\>
        Provide alternative base name for export.
        Default is <tt>meshExport</tt>.

      - \par -noBaffles
        Remove any baffles by merging the faces.

      - \par -merge
        Merge in-place interfaces

      - \par -numbered
        Use numbered patch/zone (not names) directly from ccm ids.

      - \par -remap \<name\>
        Alternative remapping dictionary
        (default: <tt>constant/remapping</tt>)

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1 (no scaling).

      - \par -solids
        Treat any solid cells present just like fluid cells.
        The default is to remove them.

Note
    - sub-domains (fluid | solid | porosity) are stored as separate domains
      within the CCM file. These are merged together to form a single mesh.
    - baffles are written as interfaces for later use

See also
    Foam::ccm::reader for more information about the File Locations

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "ccm.H"
#include "regionSplit.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reads CCM files as written by PROSTAR/STARCCM and writes an OPENFOAM "
        " polyMesh. Multi-region support for PROSTAR meshes should be stable."
        " Multi-region merging for STARCCM meshes will not always be"
        " successful."
    );

    argList::noParallel();
    argList::addArgument("ccm-file", "The input .ccm or .ccmg file");
    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of binary"
    );
    argList::addBoolOption
    (
        "export",
        "Re-export mesh in CCM format for post-processing"
    );
    argList::addBoolOption
    (
        "list",
        "List some information about the geometry"
    );
    argList::addOption
    (
        "remap",
        "name",
        "Alternative remapping dictionary (default: 'constant/remapping')"
    );
    argList::addOption
    (
        "name",
        "name",
        "Provide alternative base name when re-exporting (implies -export). "
        "Default is <meshExport>."
    );
    argList::addBoolOption
    (
        "noBaffles",
        "Remove any baffles by merging the faces"
    );
    argList::addBoolOption
    (
        "merge",
        "Merge in-place interfaces"
    );
    argList::addBoolOption
    (
        "numbered",
        "Use numbered names (eg, patch_0, zone_0) only"
    );
    argList::addOption
    (
        "scale",
        "scale",
        "Geometry scaling factor - default is 1 (ie, no scaling)"
    );
    argList::addBoolOption
    (
        "solids",
        "Treat any solid cells present just like fluid cells. "
        "The default is to remove them."
    );

    argList::noFunctionObjects();  // Never use function objects

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    runTime.functionObjects().off();  // Extra safety

    const bool optList = args.found("list");

    // exportName only has a size when export is in effect
    fileName exportName;
    if (args.readIfPresent("name", exportName))
    {
        const word ext(exportName.ext());
        // strip erroneous extension (.ccm, .ccmg, .ccmp)
        if (ext == "ccm" || ext == "ccmg" || ext == "ccmp")
        {
            exportName.remove_ext();
        }
    }
    else if (args.found("export"))
    {
        exportName = ccm::writer::defaultMeshName;
        if (args.found("case"))
        {
            exportName += '-' + args.globalCaseName();
        }
    }

    // By default, no scaling
    const scalar scaleFactor = args.getOrDefault<scalar>("scale", 1);

    // Default to binary output, unless otherwise specified
    const IOstreamOption::streamFormat format =
    (
        args.found("ascii")
      ? IOstreamOption::ASCII
      : IOstreamOption::BINARY
    );

    // More precision (for points data)
    IOstream::minPrecision(10);


    // Read control options
    // ~~~~~~~~~~~~~~~~~~~~

    ccm::reader::options rOpts;
    rOpts.removeBaffles(args.found("noBaffles"));
    rOpts.mergeInterfaces(args.found("merge"));

    if (args.found("numbered"))
    {
        rOpts.useNumberedNames(true);
    }

    if (args.found("solids"))
    {
        Info<< "treating solids like fluids" << endl;
        rOpts.keepSolid(true);
    }
    else
    {
        rOpts.keepSolid(false);
    }

    // CCM reader for reading geometry/solution
    ccm::reader reader(args.get<fileName>(1), rOpts);

    // list the geometry information
    if (optList)
    {
        Info<< "mesh geometry information:" << endl;
        if (reader.hasGeometry())
        {
            Info<< nl << "cellTable:" << reader.cellTableInfo()
                << nl << "boundaryRegion:" << reader.boundaryTableInfo()
                << nl << "interfaces:" << reader.interfaceDefinitionsInfo()
                << endl;

            if
            (
                args.found("remap")
              ? reader.remapMeshInfo(runTime, args["remap"])
              : reader.remapMeshInfo(runTime)
            )
            {
                Info<< nl
                    << "Remapped cellTable:" << reader.cellTableInfo() << nl
                    << "Remapped boundaryRegion:" << reader.boundaryTableInfo()
                    << endl;
            }
        }
        else
        {
            Info<< "NONE" << endl;
        }

        return 0;
    }
    else if (reader.readGeometry(scaleFactor))
    {
        autoPtr<polyMesh> mesh =
        (
            args.found("remap")
          ? reader.mesh(runTime, args["remap"])
          : reader.mesh(runTime)
        );

        // report mesh bounding box information
        Info<< nl << "Bounding box size: " << mesh().bounds().span() << nl;

        // check number of regions
        regionSplit rs(mesh());

        Info<< "Number of regions: " << rs.nRegions();
        if (rs.nRegions() == 1)
        {
            Info<< " (OK)." << nl;
        }
        else
        {
            Info<< nl << nl
                << "**************************************************" << nl
                << "**  WARNING: the mesh has disconnected regions  **" << nl
                << "**************************************************" << nl;
        }
        Info<< endl;
        reader.writeMesh(mesh(), format);

        // exportName only has a size when export is in effect
        if (exportName.size())
        {
            const fileName geomName = exportName + ".ccmg";
            Info<< nl << "Re-exporting geometry as " << geomName << nl;
            ccm::writer(geomName, mesh()).writeGeometry();
        }
    }
    else
    {
        FatalErrorInFunction
            << "could not read geometry"
            << exit(FatalError);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
