/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    steadyParticleTracks

Group
    grpPostProcessingUtilities

Description
    Generate a legacy VTK file of particle tracks for cases that were
    computed using a steady-state cloud

    Note:
    - Case must be re-constructed (if running in parallel) before use

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "labelPairHashes.H"
#include "IOField.H"
#include "IOobjectList.H"
#include "SortableList.H"
#include "passiveParticleCloud.H"
#include "steadyParticleTracksTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Extract list of IOobjects, modifying the input IOobjectList in the
// process
IOobjectList preFilterFields
(
    IOobjectList& cloudObjects,
    const wordRes& acceptFields,
    const wordRes& excludeFields
)
{
    IOobjectList filteredObjects(cloudObjects.capacity());
    DynamicList<label> missed(acceptFields.size());

    // Selection here is slighly different than usual
    // - an empty accept filter means "ignore everything"

    if (!acceptFields.empty())
    {
        const wordRes::filter pred(acceptFields, excludeFields);

        const wordList allNames(cloudObjects.sortedNames());

        // Detect missing fields
        forAll(acceptFields, i)
        {
            if
            (
                acceptFields[i].isLiteral()
             && !allNames.found(acceptFields[i])
            )
            {
                missed.append(i);
            }
        }

        for (const word& fldName : allNames)
        {
            const auto iter = cloudObjects.cfind(fldName);
            if (!pred(fldName) || !iter.good())
            {
                continue;  // reject
            }

            const IOobject& io = *(iter.val());

            if
            (
                //OR: fieldTypes::basic.found(io.headerClassName())
                io.isHeaderClass<IOField<label>>()
             || io.isHeaderClass<IOField<scalar>>()
             || io.isHeaderClass<IOField<vector>>()
             || io.isHeaderClass<IOField<sphericalTensor>>()
             || io.isHeaderClass<IOField<symmTensor>>()
             || io.isHeaderClass<IOField<tensor>>()
            )
            {
                // Transfer from cloudObjects -> filteredObjects
                filteredObjects.add(cloudObjects.remove(fldName));
            }
        }
    }

    if (missed.size())
    {
        WarningInFunction
            << nl
            << "Cannot find field file matching "
            << UIndirectList<wordRe>(acceptFields, missed) << endl;
    }

    return filteredObjects;
}


void readFieldsAndWriteVTK
(
    OFstream& os,
    const List<labelList>& particleMap,
    const IOobjectList& filteredObjects
)
{
    processFields<label>(os, particleMap, filteredObjects);
    processFields<scalar>(os, particleMap, filteredObjects);
    processFields<vector>(os, particleMap, filteredObjects);
    processFields<sphericalTensor>(os, particleMap, filteredObjects);
    processFields<symmTensor>(os, particleMap, filteredObjects);
    processFields<tensor>(os, particleMap, filteredObjects);
}

}  // End namespace Foam


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate a legacy VTK file of particle tracks for cases that were"
        " computed using a steady-state cloud"
    );

    argList::noParallel();
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::addOption
    (
        "dict",
        "file",
        "Alternative particleTrackDict dictionary"
    );
    argList::addVerboseOption();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const fileName vtkPath(runTime.rootPath()/runTime.globalCaseName()/"VTK");
    mkDir(vtkPath);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        const fileName vtkTimePath(vtkPath/runTime.timeName());
        mkDir(vtkTimePath);

        pointField particlePosition;
        labelList particleToTrack;
        label nTracks = 0;

        // Transfer particles to (more convenient) list
        {
            Info<< "    Reading particle positions" << endl;

            passiveParticleCloud myCloud(mesh, cloudName);
            Info<< "\n    Read " << returnReduce(myCloud.size(), sumOp<label>())
                << " particles" << endl;

            const label nParticles = myCloud.size();

            particlePosition.resize(nParticles);
            particleToTrack.resize(nParticles);

            LabelPairMap<label> trackTable;

            label np = 0;
            for (const passiveParticle& p : myCloud)
            {
                const label origId = p.origId();
                const label origProc = p.origProc();
                particlePosition[np] = p.position();

                const labelPair key(origProc, origId);

                const auto iter = trackTable.cfind(key);

                if (iter.good())
                {
                    particleToTrack[np] = *iter;
                }
                else
                {
                    particleToTrack[np] = trackTable.size();
                    trackTable.insert(key, trackTable.size());
                }

                ++np;
            }

            nTracks = trackTable.size();
        }

        if (nTracks == 0)
        {
            Info<< "\n    No track data" << endl;
        }
        else
        {
            Info<< "\n    Generating " << nTracks << " tracks" << endl;

            // Determine length of each track
            labelList trackLengths(nTracks, Zero);
            for (const label tracki : particleToTrack)
            {
                ++trackLengths[tracki];
            }

            // Particle "age" property used to sort the tracks
            List<SortableList<scalar>> agePerTrack(nTracks);
            List<List<label>> particleMap(nTracks);

            forAll(trackLengths, i)
            {
                const label length = trackLengths[i];
                agePerTrack[i].setSize(length);
                particleMap[i].setSize(length);
            }

            // Store the particle age per track
            IOobjectList cloudObjects
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudName
            );

            // TODO: gather age across all procs
            {
                tmp<IOField<scalar>> tage =
                    readParticleField<scalar>("age", cloudObjects);

                const auto& age = tage();

                labelList trackSamples(nTracks, Zero);

                forAll(particleToTrack, i)
                {
                    const label tracki = particleToTrack[i];
                    const label samplei = trackSamples[tracki];
                    agePerTrack[tracki][samplei] = age[i];
                    particleMap[tracki][samplei] = i;
                    ++trackSamples[tracki];
                }
            }


            const IOobjectList filteredObjects
            (
                preFilterFields(cloudObjects, acceptFields, excludeFields)
            );


            if (Pstream::master())
            {
                OFstream os(vtkTimePath/"particleTracks.vtk");

                Info<< "\n    Writing particle tracks to " << os.name() << endl;

                label nPoints = sum(trackLengths);

                os  << "# vtk DataFile Version 2.0" << nl
                    << "particleTracks" << nl
                    << "ASCII" << nl
                    << "DATASET POLYDATA" << nl
                    << "POINTS " << nPoints << " float" << nl;

                Info<< "\n    Writing points" << endl;

                {
                    forAll(agePerTrack, i)
                    {
                        agePerTrack[i].sort();

                        const labelList& ids = agePerTrack[i].indices();
                        labelList& particleIds = particleMap[i];

                        {
                            // Update addressing
                            List<label> sortedIds(ids);
                            forAll(sortedIds, j)
                            {
                                sortedIds[j] = particleIds[ids[j]];
                            }
                            particleIds = sortedIds;
                        }

                        forAll(ids, j)
                        {
                            const label localId = particleIds[j];
                            const point& pos = particlePosition[localId];
                            os  << pos.x() << ' ' << pos.y() << ' ' << pos.z()
                                << nl;
                        }
                    }
                }


                // Write track (line) connectivity to file

                Info<< "\n    Writing track lines" << endl;
                os  << "\nLINES " << nTracks << ' ' << nPoints + nTracks << nl;

                // Write ids of track points to file
                {
                    label globalPtI = 0;
                    forAll(particleMap, i)
                    {
                        os  << particleMap[i].size() << nl;

                        forAll(particleMap[i], j)
                        {
                            os  << ' ' << globalPtI++;

                            if (((j + 1) % 10 == 0) && (j != 0))
                            {
                                os  << nl;
                            }
                        }

                        os  << nl;
                    }
                }


                const label nFields = filteredObjects.size();

                os  << "POINT_DATA " << nPoints << nl
                    << "FIELD attributes " << nFields << nl;

                Info<< "\n    Processing fields" << nl << endl;

                readFieldsAndWriteVTK(os, particleMap, filteredObjects);
            }
        }
        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
