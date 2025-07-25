/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    particleTracks

Group
    grpPostProcessingUtilities

Description
    Generate particle tracks for cases that were computed using a
    tracked-parcel-type cloud.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "coordSetWriter.H"
#include "passiveParticleCloud.H"
#include "particleTracksSampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

template<class Type>
void writeTrackFields
(
    coordSetWriter& writer,
    HashTable<List<DynamicList<Type>>>& fieldTable
)
{
    for (const word& fieldName : fieldTable.sortedToc())
    {
        // Steal and reshape from List<DynamicList> to List<Field>
        auto& input = fieldTable[fieldName];

        List<Field<Type>> fields(input.size());
        forAll(input, tracki)
        {
            fields[tracki].transfer(input[tracki]);
        }

        writer.write(fieldName, fields);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate a file of particle tracks for cases that were"
        " computed using a tracked-parcel-type cloud"
    );
    timeSelector::addOptions();
    #include "addRegionOption.H"

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");

    argList::addOption
    (
        "dict",
        "file",
        "Alternative particleTracksProperties dictionary"
    );
    argList::addOption
    (
        "stride",
        "int",
        "Override the sample-frequency"
    );
    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to write "
        "(default: all or 'fields' from dictionary)\n"
        "Eg, 'T' or '( \"U.*\" )'"
    );
    argList::addOption
    (
        "format",
        "name",
        "The writer format "
        "(default: vtk or 'setFormat' from dictionary)"
    );
    argList::addVerboseOption();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    // ------------------------------------------------------------------------
    // Control properties

    #include "createControls.H"

    args.readListIfPresent<wordRe>("fields", acceptFields);
    args.readIfPresent("format", setFormat);

    args.readIfPresent("stride", sampleFrequency);
    sampleFrequency = max(1, sampleFrequency);  // sanity

    // Setup the writer
    auto writerPtr =
        coordSetWriter::New
        (
            setFormat,
            formatOptions.optionalSubDict(setFormat, keyType::LITERAL)
        );

    writerPtr().useTracks(true);

    if (args.verbose())
    {
        writerPtr().verbose(true);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    particleTracksSampler trackSampler;

    Info<< "Scanning times to determine track data for cloud " << cloudName
        << nl << endl;

    {
        labelList maxIds(Pstream::nProcs(), -1);
        forAll(timeDirs, timei)
        {
            runTime.setTime(timeDirs[timei], timei);
            Info<< "Time = " << runTime.timeName() << endl;

            passiveParticleCloud myCloud(mesh, cloudName);

            Info<< "    Read " << returnReduce(myCloud.size(), sumOp<label>())
                << " particles" << endl;

            for (const passiveParticle& p : myCloud)
            {
                const label origId = p.origId();
                const label origProc = p.origProc();

                // Handle case where processing particle data generated in
                // parallel using more cores than for this application.
                if (origProc >= maxIds.size())
                {
                    maxIds.resize(origProc+1, -1);
                }

                maxIds[origProc] = max(maxIds[origProc], origId);
            }
        }

        const label maxNProcs = returnReduce(maxIds.size(), maxOp<label>());
        maxIds.resize(maxNProcs, -1);

        Pstream::listCombineReduce(maxIds, maxEqOp<label>());

        // From ids to count
        const labelList numIds = maxIds + 1;

        // Set parcel addressing
        trackSampler.reset(numIds);

        Info<< nl
            << "Detected particles originating from "
            << maxNProcs << " processors." << nl
            << "Particle statistics:" << endl;

        if (Pstream::master())
        {
            const globalIndex& parcelAddr = trackSampler.parcelAddr();

            for (const label proci : parcelAddr.allProcs())
            {
                Info<< "    Found " << parcelAddr.localSize(proci)
                    << " particles originating"
                    << " from processor " << proci << nl;
            }
        }
    }

    trackSampler.setSampleRate(sampleFrequency, maxPositions, maxTracks);


    // Number of tracks to generate
    const label nTracks = trackSampler.nTracks();


    // Storage for all particle tracks
    List<DynamicList<point>> allTrackPos(nTracks);
    List<DynamicList<scalar>> allTrackTimes(nTracks);

    // Track field values by name/type
    HashTable<List<DynamicList<label>>> labelFields;  // <- mostly unused
    HashTable<List<DynamicList<scalar>>> scalarFields;
    HashTable<List<DynamicList<vector>>> vectorFields;
    HashTable<List<DynamicList<sphericalTensor>>> sphericalTensorFields;
    HashTable<List<DynamicList<symmTensor>>> symmTensorFields;
    HashTable<List<DynamicList<tensor>>> tensorFields;

    Info<< nl
        << "Generating " << nTracks
        << " particle tracks for cloud " << cloudName << nl << endl;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << " (processing)" << endl;

        // Read particles. Will be size 0 if no particles.
        passiveParticleCloud myCloud(mesh, cloudName);

        pointField localPositions(myCloud.size());
        scalarField localTimes(myCloud.size(), runTime.value());

        // Gather track data from all processors that have positions
        trackSampler.resetCloud(myCloud.size());
        {
            labelField& origIds = trackSampler.origParcelIds_;
            labelField& origProcs = trackSampler.origProcIds_;

            label np = 0;
            for (const passiveParticle& p : myCloud)
            {
                origIds[np] = p.origId();
                origProcs[np] = p.origProc();
                localPositions[np] = p.position();
                ++np;
            }

            trackSampler.gatherInplace(origIds);
            trackSampler.gatherInplace(origProcs);
        }


        // Read cloud fields (from disk) into object registry
        objectRegistry obr
        (
            IOobject
            (
                "cloudFields",
                runTime.timeName(),
                runTime
            )
        );

        myCloud.readFromFiles(obr, acceptFields, excludeFields);

        // Create track positions and track time fields
        // (not registered as IOFields)

        trackSampler.createTrackField(localPositions, allTrackPos);
        trackSampler.createTrackField(localTimes, allTrackTimes);

        // Create the track fields
        trackSampler.setTrackFields(obr, labelFields);
        trackSampler.setTrackFields(obr, scalarFields);
        trackSampler.setTrackFields(obr, vectorFields);
        trackSampler.setTrackFields(obr, sphericalTensorFields);
        trackSampler.setTrackFields(obr, symmTensorFields);
        trackSampler.setTrackFields(obr, tensorFields);
    }

    const label nFields =
    (
        labelFields.size()
      + scalarFields.size()
      + vectorFields.size()
      + sphericalTensorFields.size()
      + symmTensorFields.size()
      + tensorFields.size()
    );

    Info<< nl
        << "Extracted " << nFields << " cloud fields" << nl;

    if (nFields)
    {
        #undef  doLocalCode
        #define doLocalCode(FieldContent)                                     \
        if (!FieldContent.empty())                                            \
        {                                                                     \
            Info<< "   ";                                                     \
            for (const word& fieldName : FieldContent.sortedToc())            \
            {                                                                 \
                Info<< ' ' << fieldName;                                      \
            }                                                                 \
            Info<< nl;                                                        \
        }

        doLocalCode(labelFields);
        doLocalCode(scalarFields);
        doLocalCode(vectorFields);
        doLocalCode(sphericalTensorFields);
        doLocalCode(symmTensorFields);
        doLocalCode(tensorFields);
        #undef doLocalCode
    }

    Info<< nl
        << "Writing particle tracks (" << setFormat << " format)" << nl;

    if (Pstream::master())
    {
        PtrList<coordSet> tracks(allTrackPos.size());
        List<scalarField> times(allTrackPos.size());
        forAll(tracks, tracki)
        {
            tracks.set
            (
                tracki,
                new coordSet("track" + Foam::name(tracki), "xyz")
            );
            tracks[tracki].transfer(allTrackPos[tracki]);
            times[tracki].transfer(allTrackTimes[tracki]);

            if (!tracki) tracks[0].rename("tracks");
        }

        /// Currently unused
        /// List<vectorField> dirs(nTracks);
        /// const auto Uiter = vectorFields.cfind(UName);
        /// if (Uiter.good())
        /// {
        ///     const auto& UTracks = *Uiter;
        ///     forAll(UTracks, tracki)
        ///     {
        ///         const auto& U = UTracks[tracki];
        ///         dirs[tracki] = U/(mag(U) + ROOTVSMALL);
        ///     }
        /// }


        // Write tracks with fields
        if (nFields)
        {
            auto& writer = *writerPtr;

            const fileName outputPath
            (
                functionObject::outputPrefix/cloud::prefix/cloudName
              / "particleTracks" / "tracks"
            );

            Info<< "    "
                << runTime.relativePath(outputPath) << endl;

            writer.open(tracks, outputPath);
            writer.setTrackTimes(times);
            writer.nFields(nFields);

            writeTrackFields(writer, labelFields);
            writeTrackFields(writer, scalarFields);
            writeTrackFields(writer, vectorFields);
            writeTrackFields(writer, sphericalTensorFields);
            writeTrackFields(writer, symmTensorFields);
            writeTrackFields(writer, tensorFields);
        }
        else
        {
            Info<< "Warning: no fields, did not write" << endl;
        }
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
