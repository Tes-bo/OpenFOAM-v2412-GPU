/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    foamToVTK

Group
    grpPostProcessingUtilities

Description
    General OpenFOAM to VTK file writer.

    Other bits
    - Handles volFields, pointFields, surfaceScalarField, surfaceVectorField
      fields.
    - Mesh topo changes.
    - Output legacy or xml VTK format in ascii or binary.
    - Single time step writing.
    - Write subset only.
    - Optional decomposition of cells.

Usage
    \b foamToVTK [OPTION]

    Options:
      - \par -ascii
        Write VTK data in ASCII format instead of binary.

      - \par -legacy
        Write VTK data in legacy format instead of XML format

      - \par -fields \<fields\>
        Specify single or multiple fields to write (all by default)
        For example,
        \verbatim
          -fields T
          -fields '(p T U \"alpha.*\")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -surfaceFields
        Write surfaceScalarFields (e.g., phi)

      - \par -cellSet \<name\>
      - \par -cellZone \<name\>
        Restrict conversion to either the cellSet or the cellZone.

      - \par -faceSet \<name\>
      - \par -pointSet \<name\>
        Restrict conversion to the faceSet or pointSet.

      - \par -faceZones zone or zone list
        Specify single faceZone or or multiple faceZones (name or regex)
        to write

      - \par -nearCellValue
        Output cell value on patches instead of patch value itself

      - \par -no-boundary
        Suppress output for all boundary patches

      - \par -no-internal
        Suppress output for internal (volume) mesh

      - \par -no-lagrangian
        Suppress writing Lagrangian positions and fields.

      - \par -no-point-data
        Suppress conversion of pointFields. No interpolated PointData.

      - \par -with-ids
        Additional mesh id fields (cellID, procID, patchID)

      - \par -with-point-ids
        Additional pointID field for internal mesh

      - \par -poly-decomp
        Decompose polyhedral cells into tets/pyramids

      - \par -one-boundary
        Combine all patches into a single boundary file

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
          -exclude-patches '( inlet_1 inlet_2 "proc.*")'
        \endverbatim

Note
    The mesh subset is handled by fvMeshSubsetProxy. Slight inconsistency in
    interpolation: on the internal field it interpolates the whole volField
    to the whole-mesh pointField and then selects only those values it
    needs for the subMesh (using the fvMeshSubset cellMap(), pointMap()
    functions). For the patches however it uses the
    fvMeshSubset.interpolate function to directly interpolate the
    whole-mesh values onto the subset patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointMesh.H"
#include "emptyPolyPatch.H"
#include "volPointInterpolation.H"
#include "faceZoneMesh.H"
#include "areaFields.H"
#include "fvMeshSubsetProxy.H"
#include "faceSet.H"
#include "pointSet.H"
#include "HashOps.H"
#include "regionProperties.H"
#include "stringListOps.H"  // For stringListOps::findMatching()

#include "Cloud.H"
#include "readFields.H"
#include "reportFields.H"

#include "foamVtmWriter.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"
#include "foamVtkSurfaceMeshWriter.H"
#include "foamVtkLagrangianWriter.H"
#include "foamVtkSurfaceFieldWriter.H"
#include "foamVtkWriteTopoSet.H"
#include "foamVtkSeriesWriter.H"

#include "writeAreaFields.H"
#include "writeDimFields.H"
#include "writeVolFields.H"
#include "writePointFields.H"
#include "writeSurfaceFields.H"

#include "memInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const autoPtr<wordRes::filter>& patchSelector
)
{
    labelList indices;

    if (patchSelector && !patchSelector().empty())
    {
        // Name-based selection
        indices =
        (
            stringListOps::findMatching
            (
                patches,
                patchSelector(),
                nameOp<polyPatch>()
            )
        );
    }
    else
    {
        indices = identity(patches.size());
    }

    // Remove undesirable patches

    label count = 0;
    for (const label patchi : indices)
    {
        const polyPatch& pp = patches[patchi];

        if (isType<emptyPolyPatch>(pp))
        {
            continue;
        }
        else if (UPstream::parRun() && bool(isA<processorPolyPatch>(pp)))
        {
            break; // No processor patches for parallel output
        }

        indices[count] = patchi;
        ++count;
    }

    indices.resize(count);

    return indices;
}


//
// Process args for output options (-ascii, -legacy)
//
vtk::outputOptions getOutputOptions(const argList& args)
{
    // Default is inline ASCII xml
    vtk::outputOptions opts;

    if (args.found("legacy"))
    {
        opts.legacy(true);

        if (!args.found("ascii"))
        {
            if (sizeof(float) != 4 || sizeof(label) != 4)
            {
                opts.ascii(true);

                WarningInFunction
                    << "Using ASCII rather than legacy binary VTK format since "
                    << "float and/or label are not 4 bytes in size."
                    << nl << endl;
            }
            else
            {
                opts.ascii(false);
            }
        }
    }
    else
    {
        opts.ascii(args.found("ascii"));
    }

    return opts;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "General OpenFOAM to VTK file writer"
    );
    timeSelector::addOptions();

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");

    argList::addVerboseOption();

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of binary"
    );
    argList::addBoolOption
    (
        "legacy",
        "Write legacy format instead of xml",
        true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "poly-decomp",
        "Decompose polyhedral cells into tets/pyramids",
        true  // mark as an advanced option
    );
    argList::ignoreOptionCompat
    (
        {"xml", 1806},  // xml is now default, use -legacy to toggle
        false           // bool option, no argument
    );
    argList::ignoreOptionCompat
    (
        {"poly", 1806}, // poly is now default, use -poly-decomp to toggle
        false           // bool option, no argument
    );

    argList::addOption
    (
        "cellSet",
        "name",
        "Convert mesh subset corresponding to specified cellSet",
        true  // mark as an advanced option
    );
    argList::addOption
    (
        "cellZone",
        "name",
        "Convert mesh subset corresponding to specified cellZone",
        true  // mark as an advanced option
    );
    argList::addOption
    (
        "faceSet",
        "name",
        "Convert specified faceSet only",
        true  // mark as an advanced option
    );
    argList::addOption
    (
        "pointSet",
        "name",
        "Convert specified pointSet only",
        true  // mark as an advanced option
    );
    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple faceZones to write\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'.",
        true  // mark as an advanced option
    );
    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to write (all by default)\n"
        "Eg, 'T' or '(p T U \"alpha.*\")'"
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

    argList::addBoolOption
    (
        "processor-fields",
        "Write field values on processor boundaries only",
        true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "surfaceFields",
        "Write surfaceScalarFields (eg, phi)",
        true  // mark as an advanced option
    );
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
    argList::ignoreOptionCompat
    (
        {"finiteAreaFields", 2012},  // use -no-finite-area to disable
        false           // bool option, no argument
    );

    argList::addBoolOption
    (
        "nearCellValue",
        "Use cell value on patches instead of patch value itself",
        true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "no-boundary",
        "Suppress output for boundary patches"
    );
    argList::addBoolOption
    (
        "no-internal",
        "Suppress output for internal volume mesh"
    );
    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Suppress writing lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 1806});

    argList::addBoolOption
    (
        "no-point-data",  // noPointValues
        "Suppress conversion of pointFields. No interpolated PointData."
    );
    argList::addOptionCompat("no-point-data", {"noPointValues", 1806});

    argList::addBoolOption
    (
        "with-ids",
        "Additional mesh id fields (cellID, procID, patchID)",
        true  // mark as an advanced option
    );

    argList::addBoolOption
    (
        "with-point-ids",
        "Additional pointID field for internal mesh",
        true  // mark as an advanced option
    );

    argList::addBoolOption
    (
        "one-boundary",  // allPatches
        "Combine all patches into a single file"
    );
    argList::addOptionCompat("one-boundary", {"allPatches", 1806});

    #include "addAllRegionOptions.H"

    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to write\n"
        "Eg, 'top' or '( front \".*back\" )'"
    );
    argList::addOption
    (
        "exclude-patches",
        "wordRes",
        "Exclude single or multiple patches from writing\n"
        "Eg, 'outlet' or '( inlet \".*Wall\" )'",
        true  // mark as an advanced option
    );
    argList::addOptionCompat("exclude-patches", {"excludePatches", 2112});

    argList::ignoreOptionCompat
    (
        {"noFaceZones", 1806},  // faceZones are only enabled on demand
        false                   // bool option, no argument
    );
    argList::ignoreOptionCompat
    (
        {"noLinks", 1806},      // ignore never make any links
        false                   // bool option, no argument
    );
    argList::ignoreOptionCompat
    (
        {"useTimeName", 1806},  // ignore - works poorly with VTM formats
        false                   // bool option, no argument
    );
    argList::addBoolOption
    (
        "overwrite",
        "Remove any existing VTK output directory"
    );
    argList::addOption
    (
        "name",
        "subdir",
        "Directory name for VTK output (default: 'VTK')"
    );

    #include "setRootCase.H"

    /// const int optVerbose = args.verbose();
    const bool decomposePoly = args.found("poly-decomp");
    const bool doBoundary    = !args.found("no-boundary");
    const bool doInternal    = !args.found("no-internal");
    const bool doLagrangian  = !args.found("no-lagrangian");
    const bool doFiniteArea  = !args.found("no-finite-area");
    const bool doSurfaceFields = args.found("surfaceFields");
    const bool oneBoundary   = args.found("one-boundary") && doBoundary;
    const bool nearCellValue = args.found("nearCellValue") && doBoundary;

    const vtk::outputOptions writeOpts = getOutputOptions(args);

    bool processorFieldsOnly = false;

    if (args.found("processor-fields"))
    {
        if (!UPstream::parRun())
        {
            Info<< "Ignoring processor patch writing in serial"
                << nl << endl;
        }
        else if (writeOpts.legacy())
        {
            Info<< "Ignoring processor patch writing in legacy format"
                << nl << endl;
        }
        else
        {
            processorFieldsOnly = true;

            Info<< "Writing processor patch fields only"
                << nl << endl;
        }
    }

    if (nearCellValue)
    {
        Info<< "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }

    bool doPointValues = !args.found("no-point-data");
    if (!doPointValues)
    {
        Info<< "Point fields and interpolated point data"
            << " disabled with the '-no-point-data' option"
            << nl;
    }

    const bool withPointIds = args.found("with-point-ids");
    if (withPointIds)
    {
        Info<< "Write point ids requested";

        if (!doPointValues)
        {
            Info<< ", but ignored due to the '-no-point-data' option";
        }

        Info<< nl;
    }

    const bool withMeshIds = args.found("with-ids");
    if (withMeshIds)
    {
        Info<< "Writing mesh ids (cell, patch, proc) requested" << nl;
    }

    // Patch selection/deselection
    wordRes includedPatches, excludedPatches;
    autoPtr<wordRes::filter> patchSelector(nullptr);
    if (doBoundary)
    {
        bool resetFilter = false;
        if (args.readListIfPresent<wordRe>("patches", includedPatches))
        {
            resetFilter = true;
            Info<< "Including patches "
                << flatOutput(includedPatches) << nl << endl;
        }
        if (args.readListIfPresent<wordRe>("exclude-patches", excludedPatches))
        {
            resetFilter = true;
            Info<< "Excluding patches "
                << flatOutput(excludedPatches) << nl << endl;
        }

        if (resetFilter)
        {
            patchSelector =
                autoPtr<wordRes::filter>::New(includedPatches, excludedPatches);
        }
    }

    // Field selection/deselection
    wordRes includedFields, excludedFields;
    autoPtr<wordRes::filter> fieldSelector(nullptr);
    bool doConvertFields = !args.found("no-fields");
    if (doConvertFields)
    {
        bool resetFilter = false;
        if (args.readListIfPresent<wordRe>("fields", includedFields))
        {
            Info<< "Including fields "
                << flatOutput(includedFields) << nl << endl;

            resetFilter = !includedFields.empty();

            if (includedFields.empty())
            {
                // Compat: Can be specified as empty (ie, no fields)
                // Same as "block everything"

                doConvertFields = false;
                Info<< "Field conversion disabled by '-fields ()' option" << nl
                    << "Should use -no-fields instead" << endl;
            }
        }
        if (args.readListIfPresent<wordRe>("exclude-fields", excludedFields))
        {
            resetFilter = true;
            Info<< "Excluding fields "
                << flatOutput(excludedFields) << nl << endl;
        }

        if (resetFilter && doConvertFields)
        {
            fieldSelector =
                autoPtr<wordRes::filter>::New(includedFields, excludedFields);
        }
    }
    else if (doConvertFields)
    {
        Info<< "Field conversion disabled with the '-no-fields' option" << nl;
    }


    // Non-mandatory
    const wordRes selectedFaceZones(args.getList<wordRe>("faceZones", false));

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    // Information for file series
    HashTable<vtk::seriesWriter, fileName> vtkSeries;

    // Handle -allRegions, -regions, -region
    #include "getAllRegionOptions.H"

    // Names for sets and zones
    word cellSelectionName;
    word faceSetName;
    word pointSetName;

    fvMeshSubsetProxy::subsetType cellSubsetType = fvMeshSubsetProxy::NONE;

    string vtkName = args.globalCaseName();

    if (regionNames.size() == 1)
    {
        if (args.readIfPresent("cellSet", cellSelectionName))
        {
            vtkName = cellSelectionName;
            cellSubsetType = fvMeshSubsetProxy::SET;

            Info<< "Converting cellSet " << cellSelectionName
                << " only. New outside faces as \"oldInternalFaces\"."
                << nl;
        }
        else if (args.readIfPresent("cellZone", cellSelectionName))
        {
            vtkName = cellSelectionName;
            cellSubsetType = fvMeshSubsetProxy::ZONE;

            Info<< "Converting cellZone " << cellSelectionName
                << " only. New outside faces as \"oldInternalFaces\"."
                << nl;
        }

        args.readIfPresent("faceSet", faceSetName);
        args.readIfPresent("pointSet", pointSetName);
    }
    else
    {
        for
        (
            const char * const opt
          : { "cellSet", "cellZone", "faceSet", "pointSet" }
        )
        {
            if (args.found(opt))
            {
                Info<< "Ignoring -" << opt << " for multi-regions" << nl;
            }
        }
    }

    // ------------------------------------------------------------------------
    // Directory management

    // Sub-directory for output
    const word vtkDirName = args.getOrDefault<word>("name", "VTK");

    const fileName outputDir(args.globalPath()/vtkDirName);

    if (UPstream::master())
    {
        // Overwrite or create the VTK/regionName directories.
        // For the default region, this is simply "VTK/"

        for (const word& regionName : regionNames)
        {
            fileName regionOutDir(outputDir/polyMesh::regionName(regionName));

            if (args.found("overwrite") && Foam::isDir(regionOutDir))
            {
                Info<< "Removing old directory "
                    << args.relativePath(regionOutDir)
                    << nl << endl;
                Foam::rmDir(regionOutDir);
            }
            Foam::mkDir(regionOutDir);
        }
    }

    // ------------------------------------------------------------------------

    cpuTime timer;
    Info<< "Initial memory " << Foam::memInfo{}.size() << " kB" << endl;

    #include "createNamedMeshes.H"
    #include "createMeshAccounting.H"

    Info<< "VTK mesh topology: "
        << timer.cpuTimeIncrement() << " s, "
        << Foam::memInfo{}.size() << " kB" << endl;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        const word timeDesc = "_" + Foam::name(runTime.timeIndex());
        const scalar timeValue = runTime.value();

        Info<< "Time: " << runTime.timeName() << endl;


        // Accumulate information for multi-region VTM
        vtk::vtmWriter vtmMultiRegion;

        // vtmMultiRegion.set(vtkDir/vtkName + timeDesc)

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word& regionDir = polyMesh::regionName(regionName);

            if (regionNames.size() > 1)
            {
                Info<< "region = " << regionName << nl;
            }

            auto& meshProxy = meshProxies[regioni];
            auto& vtuMeshCells = vtuMappings[regioni];

            // polyMesh::readUpdateState meshState = mesh.readUpdate();

            // Check for new polyMesh/ and update mesh, fvMeshSubset
            // and cell decomposition.
            polyMesh::readUpdateState meshState =
                meshProxy.readUpdate();

            const fvMesh& mesh = meshProxy.mesh();

            if
            (
                meshState == polyMesh::TOPO_CHANGE
             || meshState == polyMesh::TOPO_PATCH_CHANGE
            )
            {
                // Trigger change for vtk cells too
                vtuMeshCells.clear();
            }

            // Write topoSets before attempting anything else
            {
                #include "convertTopoSet.H"
                if (wroteTopoSet)
                {
                    continue;
                }
            }

            IOobjectList objects;
            IOobjectList faObjects;

            if (doConvertFields)
            {
                // List of volume mesh objects for this instance
                objects =
                    IOobjectList(meshProxy.baseMesh(), runTime.timeName());

                // List of area mesh objects (assuming single region)
                faObjects =
                    IOobjectList
                    (
                        runTime,
                        runTime.timeName(),
                        faMesh::dbDir(meshProxy.baseMesh(), word::null),
                        IOobjectOption::NO_REGISTER
                    );

                if (fieldSelector && !fieldSelector().empty())
                {
                    objects.filterObjects(fieldSelector());
                    faObjects.filterObjects(fieldSelector());
                }

                // Remove "*_0" restart fields
                objects.prune_0();
                faObjects.prune_0();

                if (!doPointValues)
                {
                    // Prune point fields if disabled
                    objects.filterClasses
                    (
                        [](const word& clsName)
                        {
                            return fieldTypes::point.found(clsName);
                        },
                        true // prune
                    );
                }
            }

            if (processorFieldsOnly)
            {
                // Processor-patches only and continue
                #include "convertProcessorPatches.H"
                continue;
            }

            // Volume, internal, point fields
            #include "convertVolumeFields.H"

            // Surface fields
            #include "convertSurfaceFields.H"

            // Finite-area mesh and fields - need not exist
            #include "convertAreaFields.H"

            // Write lagrangian data
            #include "convertLagrangian.H"
        }

        // Emit multi-region vtm
        if (UPstream::master() && regionNames.size() > 1)
        {
            fileName outputName
            (
                outputDir/vtkName + "-regions" + timeDesc + ".vtm"
            );

            vtmMultiRegion.setTime(timeValue);
            vtmMultiRegion.write(outputName);

            fileName seriesName(vtk::seriesWriter::base(outputName));

            vtk::seriesWriter& series = vtkSeries(seriesName);

            // First time?
            // Load from file, verify against filesystem,
            // prune time >= currentTime
            if (series.empty())
            {
                series.load(seriesName, true, timeValue);
            }

            series.append(timeValue, outputName);
            series.write(seriesName);
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << Foam::memInfo{}.size() << " kB" << endl;
    }


    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << Foam::memInfo{}.peak() << " kB (peak)\n" << endl;

    return 0;
}


// ************************************************************************* //
