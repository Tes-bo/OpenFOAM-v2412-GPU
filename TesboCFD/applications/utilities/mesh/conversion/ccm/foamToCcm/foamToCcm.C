/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    foamToCcm

Group
    grpMeshConversionUtilities

Description
    Translates OPENFOAM mesh and/or results to CCM format

Usage
    \b foamToCcm [OPTION]

    Options:
      - \par -mesh
        convert mesh only to CCM format

      - \par -name \<name\>
        Provide alternative base name. Default is <tt>meshExport</tt>.

      - \par -overwrite
        No backup of existing output files.

      - \par -remap \<name\>
        Alternative remapping dictionary
        (default: <tt>constant/remapping</tt>)

      - \par -results
        Convert results only to CCM format

Note
    - No parallel data
    - No Lagrangian elements
    - the -noZero time option can be useful to avoid the often incomplete
      initial conditions (missing useful calculated values)

See also
    Foam::ccm::writer for information about the
    <tt>constant/remapping</tt> file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "volFields.H"
#include "OFstream.H"
#include "IOobjectList.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

#include "ccm.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Translate OPENFOAM data to CCM format"
    );

    timeSelector::addOptions();
    argList::noParallel();
    argList::addBoolOption
    (
        "mesh",
        "Convert mesh only"
    );
    argList::addOption
    (
        "name",
        "name",
        "Provide alternative base name. Default is <meshExport>."
    );
    argList::addBoolOption
    (
        "overwrite",
        "No backup of existing output files"
    );
    argList::addOption
    (
        "remap",
        "name",
        "Alternative remapping dictionary (default: 'constant/remapping')"
    );
    argList::addBoolOption
    (
        "results",
        "Convert results only"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"

    // The times list
    instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    const bool optMesh      = args.found("mesh");
    const bool optResults   = args.found("results");
    const bool optOverwrite = args.found("overwrite");

    fileName exportName = ccm::writer::defaultMeshName;
    if (args.readIfPresent("name", exportName))
    {
        const word ext(exportName.ext());
        // strip erroneous extension (.ccm, .ccmg, .ccmp)
        if (ext == "ccm" || ext == "ccmg" || ext == "ccmp")
        {
            exportName.remove_ext();
        }
    }
    else if (args.found("case"))
    {
        exportName += '-' + args.globalCaseName();
    }

    if (optMesh && optResults)
    {
        Warning
            << "\n-mesh and -results options are mutually exclusive\n"
            << endl;
        args.printUsage();
        FatalError.exit();
    }

//     // skip over time=0, unless some other time option has been specified
//     if
//     (
//         !args.found("noZero")
//      && !args.found("time")
//      && !args.found("latestTime")
//      && Times.size() > 2
//     )
//     {
//         startTime = 2;
//     }
//
//    runTime.setTime(Times[startTime], startTime);

    runTime.setTime(timeDirs[0], 0);
    if (optMesh)
    {
        // convert mesh only
        #include "createPolyMesh.H"

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            #include "getTimeIndex.H"

            if (timeI == 0)
            {
                ccm::writer writer
                (
                    exportName + ".ccmg",
                    mesh,
                    !optOverwrite
                );
                writer.writeGeometry();
            }
            else if (mesh.moving())
            {
                ccm::writer writer
                (
                    exportName + ".ccmg_" + timeName,
                    mesh,
                    !optOverwrite
                );
                writer.writeGeometry();
            }
        }
    }
    else
    {
        // convert fields with or without converting mesh
        #include "createNamedMesh.H"

        // #include "checkHasMovingMesh.H"
        // #include "checkHasLagrangian.H"

        IOobjectList objects(mesh, timeDirs.last().name());

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            #include "getTimeIndex.H"

            Info<< "has "
                << mesh.nCells() << " cells, "
                << mesh.nPoints() << " points, "
                << mesh.boundaryMesh().size() << " patches"
                << endl;

            if (!optResults)
            {
                if (timeI == 0)
                {
                    ccm::writer writer
                    (
                        exportName + ".ccmg",
                        mesh,
                        !optOverwrite
                    );
                    writer.writeGeometry();
                }
                else if (mesh.moving())
                {
                    ccm::writer writer
                    (
                        exportName + ".ccmg_" + timeName,
                        mesh,
                        !optOverwrite
                    );
                    writer.writeGeometry();
                }
            }

            ccm::writer writer
            (
                exportName + ".ccmp_" + timeName,
                mesh,
                !optOverwrite
            );
            // writer.setTopologyFile(exportName + ".ccmg");
            Info<< "writing solution:";
            if (args.found("remap"))
            {
                writer.writeSolution(objects, args["remap"]);
            }
            else
            {
                writer.writeSolution(objects);
            }
        }
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
