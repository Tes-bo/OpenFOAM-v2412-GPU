/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    foamToEnsight

Group
    grpPostProcessingUtilitie

Description
    Translate OpenFOAM data to EnSight format.
    An Ensight part is created for cellZones (unzoned cells are "internalMesh")
    and patches.

    - Handles volume fields, dimensioned fields, point fields
    - Handles mesh topology changes.

Usage
    \b foamToEnsight [OPTION]

    Options:
      - \par -ascii
        Write Ensight data in ASCII format instead of "C Binary"

      - \par -fields \<fields\>
        Specify single or multiple fields to write (all by default)
        For example,
        \verbatim
          -fields T
          -fields '(p T U \"alpha.*\")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -nearCellValue
        Use zero-gradient cell values on patches

      - \par -nodeValues
        Force interpolation of values to nodes

      - \par -no-boundary
        Suppress output for all boundary patches

      - \par -no-internal
        Suppress output for internal (volume) mesh

      - \par -no-cellZones
        Suppress cellZone handling

      - \par -no-lagrangian
        Suppress writing lagrangian positions and fields.

      - \par -no-mesh
        Suppress writing the geometry. Can be useful for converting partial
        results for a static geometry.

      - \par -no-point-data
        Suppress conversion of pointFields. No interpolated PointData.

      - \par -noZero
        Exclude the often incomplete initial conditions.

      - \par -index \<start\>
        Use consecutive indexing for \c data/######## files with the
        specified start index.
        Ignore the time index contained in the uniform/time file.

      - \par -name \<subdir\>
        Define sub-directory name to use for Ensight data (default: "EnSight")

      - \par -width \<n\>
        Width of Ensight data subdir (default: 8)

      - \par -cellZones NAME | LIST
        Specify single zone or multiple cell zones (name or regex) to write

      - \par -faceZones NAME | LIST
        Specify single zone or multiple face zones (name or regex) to write

      - \par -patches NAME | LIST
        Specify single patch or multiple patches (name or regex) to write
        For example,
        \verbatim
          -patches top
          -patches '( front \".*back\" )'
        \endverbatim

      - \par -exclude-patches NAME | LIST
        Exclude single or multiple patches (name or regex) from writing.
        For example,
        \verbatim
          -exclude-patches '( inlet_1 inlet_2 "proc.*" )'
        \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "Pstream.H"
#include "HashOps.H"
#include "regionProperties.H"

#include "fvc.H"
#include "fvMesh.H"
#include "fieldTypes.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "vectorIOField.H"

// file-format/conversion
#include "ensightCase.H"
#include "ensightGeoFile.H"
#include "ensightFaMesh.H"
#include "ensightMesh.H"
#include "ensightOutputCloud.H"
#include "ensightOutputAreaField.H"
#include "ensightOutputVolField.H"

// local files
#include "readFields.H"
#include "writeVolFields.H"
#include "writeDimFields.H"
#include "writePointFields.H"
#include "writeAreaFields.H"

#include "memInfo.H"

#undef foamToEnsight_useTimeIndex

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Translate OpenFOAM data to Ensight format with individual parts"
        " for cellZones, unzoned cells and patches"
    );
    timeSelector::addOptions();

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");

    argList::addVerboseOption();

    #include "addAllRegionOptions.H"

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of 'C Binary'"
    );
    argList::addOption
    (
        "index",
        "start",
        "Starting index for consecutive number of Ensight data/ files."
        " Ignore the time index contained in the uniform/time file."
        , true  // mark as an advanced option
    );
    argList::addOption
    (
        "name",
        "subdir",
        "Sub-directory name for Ensight output (default: 'EnSight')"
    );
    argList::addBoolOption
    (
        "no-overwrite",
        "Suppress removal of existing EnSight output directory"
    );
    argList::addOption
    (
        "width",
        "n",
        "Width of Ensight data subdir"
    );
    argList::addBoolOption
    (
        "nearCellValue",
        "Use zero-gradient cell values on patches"
        , true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "nodeValues",
        "Force interpolation of values to nodes"
        , true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "no-boundary",  // noPatches
        "Suppress writing any patches"
    );
    argList::addOptionCompat("no-boundary", {"noPatches", 1806});

    argList::addBoolOption
    (
        "no-internal",
        "Suppress writing the internal mesh"
    );
    argList::addBoolOption
    (
        "no-cellZones",
        "Suppress writing any cellZones"
    );
    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Suppress writing lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 1806});

    argList::addBoolOption
    (
        "no-point-data",
        "Suppress conversion of pointFields, disable -nodeValues"
    );
    argList::addBoolOption
    (
        "no-mesh", // noMesh
        "Suppress writing the geometry."
        " Can be useful for converting partial results for a static geometry"
        , true  // mark as an advanced option
    );

    // Future?
    // argList::addBoolOption
    // (
    //     "one-boundary",  // allPatches
    //     "Combine all patches into a single part"
    // );
    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress output of finite-area mesh/fields",
        true  // mark as an advanced option
    );
    argList::ignoreOptionCompat
    (
        {"finite-area", 2112},  // use -no-finite-area to disable
        false           // bool option, no argument
    );

    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to write\n"
        "Eg, 'inlet' or '(outlet \"inlet.*\")'"
    );
    argList::addOption
    (
        "exclude-patches",
        "wordRes",
        "Exclude single or multiple patches from writing\n"
        "Eg, 'outlet' or '( inlet \".*Wall\" )'"
        , true  // mark as an advanced option
    );
    argList::addOptionCompat("exclude-patches", {"excludePatches", 2112});

    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple faceZones to write\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'."
    );
    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to write (all by default)\n"
        "Eg, 'T' or '( \"U.*\" )'"
    );
    argList::addOption
    (
        "exclude-fields",
        "wordRes",
        "Exclude single or multiple fields",
        true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "no-fields",
        "Suppress conversion of fields"
    );

    argList::addOption
    (
        "cellZones",
        "wordRes",
        "Specify single or multiple cellZones to write\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'."
    );
    argList::addOptionCompat("cellZones", {"cellZone", 1912});

    #include "setRootCase.H"

    // ------------------------------------------------------------------------
    // Configuration

    // Default to binary output, unless otherwise specified
    const IOstreamOption::streamFormat format =
    (
        args.found("ascii")
      ? IOstreamOption::ASCII
      : IOstreamOption::BINARY
    );

    const int optVerbose = args.verbose();
    const bool doBoundary    = !args.found("no-boundary");
    const bool doInternal    = !args.found("no-internal");
    const bool doCellZones   = !args.found("no-cellZones");
    const bool doLagrangian  = !args.found("no-lagrangian");
    const bool doFiniteArea  = !args.found("no-finite-area");
    const bool doPointValues = !args.found("no-point-data");
    const bool nearCellValue = args.found("nearCellValue") && doBoundary;

    // Control for numbering iterations
    label indexingNumber(0);
    const bool doConsecutive = args.readIfPresent("index", indexingNumber);

    // Write the geometry, unless otherwise specified
    bool doGeometry = !args.found("no-mesh");

    if (nearCellValue)
    {
        Info<< "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }
    if (!doPointValues)
    {
        Info<< "Point fields and interpolated point data"
            << " disabled with the '-no-point-data' option"
            << nl;
    }

    //
    // General (case) output options
    //
    ensightCase::options caseOpts(format);

    // Forced point interpolation?
    caseOpts.nodeValues(doPointValues && args.found("nodeValues"));
    caseOpts.width(args.getOrDefault<label>("width", 8));
    caseOpts.overwrite(!args.found("no-overwrite")); // Remove existing?

    // Can also have separate directory for lagrangian
    // caseOpts.separateCloud(true);

    ensightMesh::options writeOpts;
    writeOpts.useBoundaryMesh(doBoundary);
    writeOpts.useInternalMesh(doInternal);
    writeOpts.useCellZones(doCellZones);

    // Patch selection/deselection
    if (args.found("patches"))
    {
        writeOpts.patchSelection(args.getList<wordRe>("patches"));
    }
    if (args.found("exclude-patches"))
    {
        writeOpts.patchExclude(args.getList<wordRe>("exclude-patches"));
    }

    if (args.found("faceZones"))
    {
        writeOpts.faceZoneSelection(args.getList<wordRe>("faceZones"));
    }
    if (args.found("cellZones"))
    {
        writeOpts.cellZoneSelection(args.getList<wordRe>("cellZones"));
    }

    // Report the setup
    writeOpts.print(Info);

    // Field selection/deselection
    wordRes includedFields, excludedFields;
    autoPtr<wordRes::filter> fieldSelector(nullptr);
    const bool doConvertFields = !args.found("no-fields");
    if (doConvertFields)
    {
        bool resetFilter = false;
        if (args.readListIfPresent<wordRe>("fields", includedFields))
        {
            resetFilter = true;
            Info<< "Including fields "
                << flatOutput(includedFields) << nl << endl;
        }
        if (args.readListIfPresent<wordRe>("exclude-fields", excludedFields))
        {
            resetFilter = true;
            Info<< "Excluding fields "
                << flatOutput(excludedFields) << nl << endl;
        }
        if (resetFilter)
        {
            fieldSelector =
                autoPtr<wordRes::filter>::New(includedFields, excludedFields);
        }
    }
    else if (doConvertFields)
    {
        Info<< "Field conversion disabled with the '-no-fields' option" << nl;
    }

    // ------------------------------------------------------------------------

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    // Handle -allRegions, -regions, -region
    #include "getAllRegionOptions.H"

    // ------------------------------------------------------------------------
    // Directory management

    // Define sub-directory name to use for EnSight data.
    // The path to the ensight directory is at case level only
    // - For parallel cases, data only written from master

    // Sub-directory for output
    const word ensDirName = args.getOrDefault<word>("name", "EnSight");

    fileName outputDir(args.globalPath()/ensDirName);

    if (!outputDir.isAbsolute())
    {
        outputDir = args.globalPath()/outputDir;
    }


    // ------------------------------------------------------------------------
    cpuTime timer;
    Info<< "Initial memory " << Foam::memInfo{}.size() << " kB" << endl;

    #include "createNamedMeshes.H"
    #include "createMeshAccounting.H"

    if (Pstream::master())
    {
        Info<< "Converting " << timeDirs.size() << " time steps" << nl;
        // ensCase.printInfo(Info) << endl;
    }

    // Check mesh motion
    #include "checkMeshMoving.H"
    if (hasMovingMesh && !doGeometry)
    {
        Info<< "has moving mesh: ignoring '-no-mesh' option" << endl;
        doGeometry = true;
    }

    // Check lagrangian
    #include "findCloudFields.H"

    // Check field availability
    #include "checkFieldAvailability.H"

    // test the pre-check variable if there is a moving mesh
    // time-set for geometries
    // TODO: split off into separate time-set,
    // but need to verify ensight spec

    Info<< "Startup in "
        << timer.cpuTimeIncrement() << " s, "
        << Foam::memInfo{}.size() << " kB" << nl << endl;


    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        // Index for the Ensight case(s). Continues if not possible
        #include "getTimeIndex.H"

        Info<< "Time [" << timeIndex << "] = " << runTime.timeName() << nl;

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word& regionDir = polyMesh::regionName(regionName);

            if (regionNames.size() > 1)
            {
                Info<< "region=" << regionName << nl;
            }

            auto& mesh = meshes[regioni];

            polyMesh::readUpdateState meshState = mesh.readUpdate();
            const bool moving = (meshState != polyMesh::UNCHANGED);

            // Ensight
            auto& ensCase = ensightCases[regioni];
            auto& ensMesh = ensightMeshes[regioni];

            // Finite-area (can be missing)
            auto* ensFaCasePtr = ensightCasesFa.get(regioni);
            auto* ensFaMeshPtr = ensightMeshesFa.get(regioni);

            ensCase.setTime(timeDirs[timei], timeIndex);
            if (ensFaCasePtr)
            {
                ensFaCasePtr->setTime(timeDirs[timei], timeIndex);
            }

            if (moving)
            {
                ensMesh.expire();
                ensMesh.correct();

                if (ensFaMeshPtr)
                {
                    ensFaMeshPtr->expire();
                    ensFaMeshPtr->correct();
                }
            }

            if ((timei == 0 || moving) && doGeometry)
            {
                // finite-volume
                {
                    autoPtr<ensightGeoFile> os =
                        ensCase.newGeometry(hasMovingMesh);

                    ensMesh.write(os.ref());
                }

                // finite-area
                if (ensFaCasePtr && ensFaMeshPtr)
                {
                    autoPtr<ensightGeoFile> os =
                        ensFaCasePtr->newGeometry(hasMovingMesh);

                    ensFaMeshPtr->write(os.ref());
                }
            }

            // Objects at this time
            IOobjectList objects(mesh, runTime.timeName());

            objects.filterObjects
            (
                availableRegionObjectNames[regioni]
            );

            // Volume, internal, point fields
            #include "convertVolumeFields.H"

            // The finite-area objects at this time
            IOobjectList faObjects;

            if (ensFaMeshPtr)
            {
                faObjects =
                    IOobjectList(ensFaMeshPtr->mesh(), runTime.timeName());

                faObjects.filterObjects
                (
                    availableFaRegionObjectNames[regioni]
                );
            }

            // The finiteArea fields
            #include "convertAreaFields.H"

            // Lagrangian fields
            #include "convertLagrangian.H"
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << Foam::memInfo{}.size() << " kB" << nl << nl;
    }

    // Write cases
    forAll(ensightCases, regioni)
    {
        ensightCases[regioni].write();
    }

    forAll(ensightCasesFa, regioni)
    {
        if (ensightCasesFa.set(regioni))
        {
            ensightCasesFa[regioni].write();
        }
    }

    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << Foam::memInfo{}.peak() << " kB (peak)" << nl << endl;

    return 0;
}


// ************************************************************************* //
