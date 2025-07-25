/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
    checkMesh

Group
    grpMeshManipulationUtilities

Description
    Checks validity of a mesh.

Usage
    \b checkMesh [OPTION]

    Options:
      - \par -allGeometry
        Checks all (including non finite-volume specific) geometry

      - \par -allTopology
        Checks all (including non finite-volume specific) addressing

      - \par -meshQuality
        Checks against user defined (in \a system/meshQualityDict) quality
        settings

      - \par -region \<name\>
        Specify an alternative mesh region.

      - \par -allRegions
        Check all regions in regionProperties.

    \param -writeSets \<surfaceFormat\> \n
    Reconstruct all cellSets and faceSets geometry and write to postProcessing
    directory according to surfaceFormat (e.g. vtk or ensight). Additionally
    reconstructs all pointSets and writes as vtk format.

    \param -writeAllFields \n
    Writes all mesh quality measures as fields.

    \param -writeAllSurfaceFields \n
    Adds writing of surface fields when used in combination with writeAllFields.

    \param -writeFields '(\<fieldName\>)' \n
    Writes selected mesh quality measures as fields.

    \param -write-edges \n
    Write bad edges (possibly relevant for finite-area) in VTK format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "globalMeshData.H"
#include "vtkCoordSetWriter.H"
#include "vtkSurfaceWriter.H"
#include "IOdictionary.H"
#include "regionProperties.H"
#include "polyMeshTools.H"

#include "checkTools.H"
#include "checkTopology.H"
#include "checkGeometry.H"
#include "checkMeshQuality.H"
#include "writeFields.H"

#include "OFstream.H"
#include "JSONformatter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Checks validity of a mesh"
    );

    timeSelector::addOptions();
    #include "addAllRegionOptions.H"

    argList::addBoolOption
    (
        "noTopology",
        "Skip checking the mesh topology"
    );
    argList::addBoolOption
    (
        "allGeometry",
        "Include bounding box checks"
    );
    argList::addBoolOption
    (
        "allTopology",
        "Include extra topology checks"
    );
    argList::addBoolOption
    (
        "writeAllFields",
        "Write volFields with mesh quality parameters"
    );
    argList::addBoolOption
    (
        "writeAllSurfaceFields",
        "Write surfaceFields with mesh quality parameters"
    );
    argList::addOption
    (
        "writeFields",
        "wordList",
        "Write volFields with selected mesh quality parameters"
    );
    argList::addBoolOption
    (
        "meshQuality",
        "Read user-defined mesh quality criteria from system/meshQualityDict"
    );
    argList::addOption
    (
        "writeSets",
        "surfaceFormat",
        "Reconstruct and write all faceSets and cellSets in selected format"
    );
    argList::addBoolOption
    (
        "write-edges",
        "Write bad edges (possibly relevant for finite-area) in vtk format"
    );

    argList::addOption
    (
        "writeChecks",
        "word",
        "Write checks to file in dictionary or JSON format"
    );


    #include "setRootCase.H"
    #include "createTime.H"
    #include "getAllRegionOptions.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMeshes.H"

    const bool noTopology  = args.found("noTopology");
    const bool allGeometry = args.found("allGeometry");
    const bool allTopology = args.found("allTopology");
    const bool meshQuality = args.found("meshQuality");
    const bool optWriteEdges = args.found("write-edges");

    const word surfaceFormat = args.getOrDefault<word>("writeSets", "");
    const bool writeSets = surfaceFormat.size();


    // All potential writeable fields
    const wordHashSet allFields
    ({
        "nonOrthoAngle",
        "faceWeight",
        "skewness",
        "cellDeterminant",
        "aspectRatio",
        "cellShapes",
        "cellVolume",
        "cellVolumeRatio",
        "cellAspectRatio",
        "minTetVolume",
        "minPyrVolume",
        "cellRegion",
        "wallDistance",
        "cellZone",
        "faceZone"
    });

    #include "writeMeshChecks.H"

    if (args.found("writeChecks"))
    {
        writeChecksFormat =
            writeChecksFormatTypeNames.get(args.get<word>("writeChecks"));
    }

    const bool writeFaceFields = args.found("writeAllSurfaceFields");
    wordHashSet selectedFields;
    if (args.found("writeFields"))
    {
        selectedFields = args.getList<word>("writeFields");
        const wordHashSet badFields(selectedFields - allFields);

        if (!badFields.empty())
        {
            FatalErrorInFunction
                << "Illegal field(s): "
                << flatOutput(badFields.sortedToc()) << nl
                << "Valid fields: "
                << flatOutput(allFields.sortedToc()) << nl
                << exit(FatalError);
        }
    }
    else if (args.found("writeAllFields"))
    {
        selectedFields = allFields;
    }
    else if (writeFaceFields)
    {
        FatalErrorInFunction
            << "Option 'writeAllSurfaceFields' only valid in combination"
            << " with 'writeFields' or 'writeAllFields'"
            << nl << exit(FatalError);
    }


    Info<< "Check mesh..." << nl;
    Info().incrIndent();

    if (noTopology)
    {
        Info<< indent
            << "Disabling all topology checks." << nl;
    }
    if (allTopology)
    {
        Info<< indent
            << "Enabling all (cell, face, edge, point) topology checks." << nl;
    }
    if (allGeometry)
    {
        Info<< indent
            << "Enabling all geometry checks." << nl;
    }
    if (meshQuality)
    {
        Info<< indent
            << "Enabling user-defined geometry checks." << nl;
    }
    if (writeSets)
    {
        Info<< indent
            << "Reconstructing and writing " << surfaceFormat
            << " representation of all faceSets and cellSets." << nl;
    }
    if (selectedFields.size())
    {
        Info<< indent
            << "Writing mesh quality as fields "
            << flatOutput(selectedFields.sortedToc()) << nl;
    }
    if (optWriteEdges)
    {
        Info<< indent
            << "Writing any bad edges in vtk format" << nl;
    }
    Info().decrIndent();
    Info<< endl;


    PtrList<IOdictionary> qualDict(meshes.size());
    if (meshQuality)
    {
        forAll(meshes, meshi)
        {
            qualDict.set
            (
                meshi,
                new IOdictionary
                (
                    IOobject
                    (
                        "meshQualityDict",
                        meshes[meshi].time().system(),
                        meshes[meshi],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );
        }
    }


    autoPtr<surfaceWriter> surfWriter;
    autoPtr<coordSetWriter> setWriter;
    if (writeSets)
    {
        surfWriter = surfaceWriter::New(surfaceFormat);
        setWriter = coordSetWriter::New(coordSetWriters::vtkWriter::typeName);
    }


    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        // Get most changed of all meshes
        polyMesh::readUpdateState state = polyMesh::UNCHANGED;
        for (auto& mesh : meshes)
        {
            state = polyMeshTools::combine(state, mesh.readUpdate());
        }


        if
        (
            !timei
         || state == polyMesh::TOPO_CHANGE
         || state == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            forAll(meshes, meshi)
            {
                const auto& mesh = meshes[meshi];

                // Reconstruct globalMeshData
                mesh.globalData();

                printMeshStats(mesh, allTopology);

                label nFailedChecks = 0;

                if (!noTopology)
                {
                    nFailedChecks += checkTopology
                    (
                        mesh,
                        allTopology,
                        allGeometry,
                        surfWriter,
                        setWriter,
                        optWriteEdges
                    );
                }

                nFailedChecks += checkGeometry
                (
                    mesh,
                    allGeometry,
                    surfWriter,
                    setWriter
                );

                if (meshQuality)
                {
                    nFailedChecks +=
                        checkMeshQuality(mesh, qualDict[meshi], surfWriter);
                }


                // Note: no reduction in nFailedChecks necessary since is
                //       counter of checks, not counter of failed cells,faces
                //       etc.

                if (nFailedChecks == 0)
                {
                    Info<< "\nMesh OK.\n" << endl;
                }
                else
                {
                    Info<< "\nFailed " << nFailedChecks << " mesh checks.\n"
                        << endl;
                }


                // Write selected fields
                Foam::writeFields(mesh, selectedFields, writeFaceFields);

                writeMeshChecks(mesh, writeChecksFormat);
            }
        }
        else if (state == polyMesh::POINTS_MOVED)
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            forAll(meshes, meshi)
            {
                const auto& mesh = meshes[meshi];

                label nFailedChecks = checkGeometry
                (
                    mesh,
                    allGeometry,
                    surfWriter,
                    setWriter
                );

                if (meshQuality)
                {
                    nFailedChecks +=
                        checkMeshQuality(mesh, qualDict[meshi], surfWriter);
                }


                if (nFailedChecks)
                {
                    Info<< "\nFailed " << nFailedChecks << " mesh checks.\n"
                        << endl;
                }
                else
                {
                    Info<< "\nMesh OK.\n" << endl;
                }


                // Write selected fields
                Foam::writeFields(mesh, selectedFields, writeFaceFields);

                writeMeshChecks(mesh, writeChecksFormat);
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
