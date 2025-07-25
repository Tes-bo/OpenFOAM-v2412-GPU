/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "surfaceNoise.H"
#include "surfaceReader.H"
#include "surfaceWriter.H"
#include "globalIndex.H"
#include "argList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace noiseModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceNoise, 0);
addToRunTimeSelectionTable(noiseModel, surfaceNoise, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void surfaceNoise::initialise(const fileName& fName)
{
    Info<< "Reading data file: "
        << fileObr_.time().relativePath(fName) << endl;

    instantList allTimes;
    label nAvailableTimes = 0;

    // All reading performed on the master processor only
    if (Pstream::master())
    {
        // Create the surface reader
        readerPtr_ = surfaceReader::New(readerType_, fName);

        // Find the index of the pressure data
        const wordList fieldNames(readerPtr_->fieldNames(0));
        pIndex_ = fieldNames.find(pName_);
        if (pIndex_ == -1)
        {
            FatalErrorInFunction
                << "Unable to find pressure field name " << pName_
                << " in list of available fields: " << fieldNames
                << exit(FatalError);
        }

        // Create the surface writer
        // - Could be done later, but since this utility can process a lot of
        //   data we can ensure that the user-input is correct prior to doing
        //   the heavy lifting

        // Set the time range
        allTimes = readerPtr_->times();
        startTimeIndex_ = instant::findStart(allTimes, startTime_);

        // Determine the windowing
        nAvailableTimes = allTimes.size() - startTimeIndex_;
    }

    Pstream::broadcasts
    (
        UPstream::worldComm,
        pIndex_,
        startTimeIndex_,
        nAvailableTimes
    );


    // Note: all processors should call the windowing validate function
    label nRequiredTimes = windowModelPtr_->validate(nAvailableTimes);

    if (UPstream::master())
    {
        // Restrict times
        times_.resize_nocopy(nRequiredTimes);
        forAll(times_, timei)
        {
            times_[timei] = allTimes[timei + startTimeIndex_].value();
        }

        deltaT_ =
        (
            sampleFreq_ > 0
          ? (1.0/sampleFreq_)
          : checkUniformTimeStep(times_)
        );

        // Read the surface geometry
        // Note: hard-coded to read mesh from first time index
        const meshedSurface& surf = readerPtr_->geometry(0);

        nFaces_ = surf.nFaces();
    }

    Pstream::broadcasts
    (
        UPstream::worldComm,
        times_,
        deltaT_,
        nFaces_
    );
}


void surfaceNoise::readSurfaceData
(
    const globalIndex& procFaceAddr,
    List<scalarField>& pData
)
{
    // Data is stored as pressure on surface at a given time.  Now we need to
    // pivot the data so that we have the complete pressure time history per
    // surface face.  In serial mode, this results in all pressure data being
    // loaded into memory (!)

    const label nLocalFace = procFaceAddr.localSize();

    // Complete pressure time history data for subset of faces
    pData.resize_nocopy(nLocalFace);
    const label nTimes = times_.size();
    for (scalarField& pf : pData)
    {
        pf.resize_nocopy(nTimes);
    }

    Info<< "Reading pressure data" << endl;

    // Master only
    scalarField allData;

    if (Pstream::parRun())
    {
        // Procedure:
        // 1. Master processor reads pressure data for all faces for all times
        // 2. Master sends each processor a subset of faces
        // 3. Local processor reconstructs the full time history for the subset
        //    of faces
        // Note: reading all data on master to avoid potential NFS problems...

        scalarField scratch;

        if (!useBroadcast_)
        {
            scratch.resize(nLocalFace);
        }

        // Read data and send to sub-ranks
        forAll(times_, timei)
        {
            const label fileTimeIndex = timei + startTimeIndex_;

            if (Pstream::master())
            {
                Info<< "    time: " << times_[timei] << endl;

                // Read pressure at all faces for time timeI
                allData = readerPtr_->field(fileTimeIndex, pIndex_, scalar(0));
            }

            if (useBroadcast_)
            {
                Pstream::broadcast(allData);
            }
            else
            {
                procFaceAddr.scatter
                (
                    allData,
                    scratch,
                    UPstream::msgType(),
                    commType_,
                    UPstream::worldComm
                );
            }

            scalarField::subField procData =
            (
                useBroadcast_
              ? allData.slice(procFaceAddr.range())
              : scratch.slice(0, nLocalFace)
            );

            // Apply conversions
            procData *= rhoRef_;

            // Transcribe this time snapshot (transpose)
            forAll(procData, facei)
            {
                pData[facei][timei] = procData[facei];
            }
        }
    }
    else
    {
        // Read data - no sub-ranks
        forAll(times_, timei)
        {
            const label fileTimeIndex = timei + startTimeIndex_;

            Info<< "    time: " << times_[timei] << endl;

            allData = readerPtr_->field(fileTimeIndex, pIndex_, scalar(0));

            auto& procData = allData;

            // Apply conversions
            procData *= rhoRef_;

            // Transcribe this time snapshot (transpose)
            forAll(procData, facei)
            {
                pData[facei][timei] = procData[facei];
            }
        }
    }

    forAll(pData, facei)
    {
        pData[facei] -= average(pData[facei]);
    }


    Info<< "Read "
        << returnReduce(pData.size(), sumOp<label>())
        << " pressure traces with "
        << times_.size()
        << " time values" << nl << endl;
}


scalar surfaceNoise::surfaceAverage
(
    const scalarField& data,
    const globalIndex& procFaceAddr
) const
{
    if (!nFaces_)
    {
        // Already reduced, can use as sanity check
        return 0;
    }

    scalar areaAverage = 0;

    if (areaAverage_)
    {
        if (Pstream::parRun())
        {
            // Collect the surface data so that we can output the surfaces
            scalarField allData;

            procFaceAddr.gather
            (
                data,
                allData,
                UPstream::msgType(),
                commType_,
                UPstream::worldComm
            );

            if (Pstream::master())
            {
                // Note: hard-coded to read mesh from first time index
                const meshedSurface& surf = readerPtr_->geometry(0);

                areaAverage = sum(allData*surf.magSf())/sum(surf.magSf());
            }
        }
        else
        {
            // Note: hard-coded to read mesh from first time index
            const meshedSurface& surf = readerPtr_->geometry(0);

            areaAverage = sum(data*surf.magSf())/sum(surf.magSf());
        }

        Pstream::broadcast(areaAverage);
    }
    else
    {
        // Ensemble averaged
        // - same as gAverage, but already know number of faces

        areaAverage = sum(data);
        reduce(areaAverage, sumOp<scalar>());

        areaAverage /= (scalar(nFaces_) + ROOTVSMALL);
    }

    return areaAverage;
}


scalar surfaceNoise::writeSurfaceData
(
    const fileName& outDirBase,
    const word& fName,
    const word& title,
    const scalar freq,
    const scalarField& data,
    const globalIndex& procFaceAddr,
    const bool writeSurface
) const
{
    Info<< "    processing " << title << " for frequency " << freq << endl;

    const instant freqInst(freq, Foam::name(freq));

    if (!writeSurface)
    {
        return surfaceAverage(data, procFaceAddr);
    }

    scalar areaAverage = 0;

    if (Pstream::parRun())
    {
        // Collect the surface data so that we can output the surfaces
        scalarField allData;

        procFaceAddr.gather
        (
            data,
            allData,
            UPstream::msgType(),
            commType_,
            UPstream::worldComm
        );

        if (Pstream::master())
        {
            // Note: hard-coded to read mesh from first time index
            const meshedSurface& surf = readerPtr_->geometry(0);

            if (areaAverage_)
            {
                areaAverage = sum(allData*surf.magSf())/sum(surf.magSf());
            }
            else
            {
                areaAverage = sum(allData)/(allData.size() + ROOTVSMALL);
            }

            // (writeSurface == true)
            {
                // Time-aware, with time spliced into the output path
                writerPtr_->beginTime(freqInst);

                writerPtr_->open
                (
                    surf.points(),
                    surf.surfFaces(),
                    (outDirBase / fName),
                    false  // serial - already merged
                );

                writerPtr_->nFields(1); // Legacy VTK
                writerPtr_->write(title, allData);

                writerPtr_->endTime();
                writerPtr_->clear();
            }
        }
    }
    else
    {
        // Note: hard-coded to read mesh from first time index
        const meshedSurface& surf = readerPtr_->geometry(0);

        if (areaAverage_)
        {
            areaAverage = sum(data*surf.magSf())/sum(surf.magSf());
        }
        else
        {
            areaAverage = sum(data)/(data.size() + ROOTVSMALL);
        }

        // (writeSurface == true)
        {
            // Time-aware, with time spliced into the output path
            writerPtr_->beginTime(freqInst);

            writerPtr_->open
            (
                surf.points(),
                surf.surfFaces(),
                (outDirBase / fName),
                false  // serial - already merged
            );

            writerPtr_->nFields(1); // Legacy VTK
            writerPtr_->write(title, data);

            writerPtr_->endTime();
            writerPtr_->clear();
        }
    }

    Pstream::broadcast(areaAverage);
    return areaAverage;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceNoise::surfaceNoise
(
    const dictionary& dict,
    const objectRegistry& obr,
    const word& name,
    const bool readFields
)
:
    noiseModel(dict, obr, name, false),
    inputFileNames_(),
    pName_("p"),
    pIndex_(0),
    times_(),
    deltaT_(0),
    startTimeIndex_(0),
    nFaces_(0),
    fftWriteInterval_(1),
    areaAverage_(false),
    useBroadcast_(false),
    commType_(UPstream::commsTypes::scheduled),
    readerType_(),
    readerPtr_(nullptr),
    writerPtr_(nullptr)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool surfaceNoise::read(const dictionary& dict)
{
    if (noiseModel::read(dict))
    {
        if (!dict.readIfPresent("files", inputFileNames_))
        {
            inputFileNames_.resize(1);
            dict.readEntry("file", inputFileNames_.first());
        }

        dict.readIfPresent("p", pName_);
        dict.readIfPresent("fftWriteInterval", fftWriteInterval_);

        Info<< this->type() << nl
            << "    Pressure field name: " << pName_ << nl
            << "    FFT write interval: " << fftWriteInterval_ << nl;

        dict.readIfPresent("areaAverage", areaAverage_);

        if (areaAverage_)
        {
            Info<< "    Averaging: area weighted" << endl;
        }
        else
        {
            Info<< "    Averaging: ensemble" << endl;
        }

        useBroadcast_ = false;
        dict.readIfPresent("broadcast", useBroadcast_);
        UPstream::commsTypeNames.readIfPresent("commsType", dict, commType_);

        if (Pstream::parRun())
        {
            Info<< "    Distribute fields: "
                << UPstream::commsTypeNames[commType_];

            if (useBroadcast_)
            {
                Info<< " (broadcast all)";
            }
            Info<< endl;
        }

        readerType_ = dict.get<word>("reader");

        // Surface writer (keywords: writer, writeOptions)

        const word writerType = dict.get<word>("writer");

        writerPtr_ = surfaceWriter::New
        (
            writerType,
            surfaceWriter::formatOptions(dict, writerType, "writeOptions")
        );

        // Use outputDir/TIME/surface-name
        writerPtr_->useTimeDir(true);

        Info<< endl;

        return true;
    }

    return false;
}


void surfaceNoise::calculate()
{
    forAll(inputFileNames_, filei)
    {
        fileName fName = inputFileNames_[filei];
        fName.expand();

        if (!fName.isAbsolute())
        {
            fName = argList::envGlobalPath()/fName;
        }

        initialise(fName);

        // Processor face addressing
        globalIndex procFaceAddr;

        if (Pstream::parRun())
        {
            // Calculate face/proc offsets manually
            labelList procFaceOffsets(Pstream::nProcs() + 1);
            const label nFacePerProc = floor(nFaces_/Pstream::nProcs()) + 1;

            procFaceOffsets[0] = 0;
            for (label proci = 1; proci < procFaceOffsets.size(); ++proci)
            {
                procFaceOffsets[proci] = min(proci*nFacePerProc, nFaces_);
            }

            procFaceAddr.offsets() = std::move(procFaceOffsets);

            // Don't need to broadcast. Already identical on all ranks
        }
        else
        {
            // Local data only
            procFaceAddr.reset(globalIndex::gatherNone{}, nFaces_);
        }

        // Pressure time history data per face
        List<scalarField> pData;

        // Read pressure data from file
        readSurfaceData(procFaceAddr, pData);

        // Process the pressure data, and store results as surface values per
        // frequency so that it can be output using the surface writer

        Info<< "Creating noise FFTs" << endl;

        const scalarField freq1(uniformFrequencies(deltaT_, true));

        const scalar maxFreq1 = max(freq1);

        // Reset desired frequency range if outside actual frequency range
        fLower_ = min(fLower_, maxFreq1);
        fUpper_ = min(fUpper_, maxFreq1);

        // Storage for FFT data
        const label nLocalFace = pData.size();
        const label nFFT = ceil(freq1.size()/scalar(fftWriteInterval_));

        List<scalarField> surfPrmsf(nFFT);
        List<scalarField> surfPSDf(nFFT);
        forAll(surfPrmsf, freqI)
        {
            surfPrmsf[freqI].setSize(nLocalFace);
            surfPSDf[freqI].setSize(nLocalFace);
        }

        // Storage for 1/3 octave data
        labelList octave13BandIDs;
        scalarField octave13FreqCentre;
        setOctaveBands
        (
            freq1,
            fLower_,
            fUpper_,
            3,
            octave13BandIDs,
            octave13FreqCentre
        );

        label bandSize = 0;
        if (octave13BandIDs.empty())
        {
            WarningInFunction
                << "Octave band calculation failed (zero sized). "
                << "Please check your input data"
                << endl;
        }
        else
        {
            bandSize = octave13BandIDs.size() - 1;
        }

        List<scalarField> surfPrms13f(bandSize);
        forAll(surfPrms13f, freqI)
        {
            surfPrms13f[freqI].setSize(nLocalFace);
        }

        const windowModel& win = windowModelPtr_();

        {
            forAll(pData, faceI)
            {
                const scalarField& p = pData[faceI];

                // Generate the FFT-based data
                const scalarField Prmsf(RMSmeanPf(p));
                const scalarField PSDf(this->PSDf(p, deltaT_));

                // Store the frequency results in slot for face of surface
                forAll(surfPrmsf, i)
                {
                    label freqI = i*fftWriteInterval_;
                    surfPrmsf[i][faceI] = Prmsf[freqI];
                    surfPSDf[i][faceI] = PSDf[freqI];
                }

                // Integrated PSD = P(rms)^2 [Pa^2]
                const scalarField Prms13f
                (
                    octaves
                    (
                        PSDf,
                        freq1,
                        octave13BandIDs
                    )
                );

                // Store the 1/3 octave results in slot for face of surface
                forAll(surfPrms13f, freqI)
                {
                    surfPrms13f[freqI][faceI] = Prms13f[freqI];
                }
            }
        }

        const word fNameBase = fName.stem();

        // Output directory
        fileName outDirBase(baseFileDir(filei)/fNameBase);

        const scalar deltaf = 1.0/(deltaT_*win.nSamples());
        Info<< "Writing fft surface data";
        if (fftWriteInterval_ == 1)
        {
            Info<< endl;
        }
        else
        {
            Info<< " at every " << fftWriteInterval_ << " frequency points"
                << endl;
        }

        // Common output information
        // Note: hard-coded to read mesh from first time index
        scalar surfArea = 0;
        label surfSize = 0;
        if (Pstream::master())
        {
            const meshedSurface& surf = readerPtr_->geometry(0);
            surfArea = sum(surf.magSf());
            surfSize = surf.size();
        }
        Pstream::broadcasts
        (
            UPstream::worldComm,
            surfArea,
            surfSize
        );

        List<Tuple2<string, token>> commonInfo
        ({
            {"Area average", token(word(Switch::name(areaAverage_)))},
            {"Area sum", token(surfArea)},
            {"Number of faces", token(surfSize)}
        });

        {
            fileName outDir(outDirBase/"fft");
            fileName outSurfDir(filePath(outDir));

            // Determine frequency range of interest
            // Note: frequencies have fixed interval, and are in the range
            //       0 to fftWriteInterval_*(n-1)*deltaf
            label f0 = ceil(fLower_/deltaf/scalar(fftWriteInterval_));
            label f1 = floor(fUpper_/deltaf/scalar(fftWriteInterval_));
            label nFreq = f1 - f0;

            scalarField PrmsfAve(nFreq, Zero);
            scalarField PSDfAve(nFreq, Zero);
            scalarField fOut(nFreq, Zero);

            if (nFreq == 0)
            {
                WarningInFunction
                    << "No surface data available using a fftWriteInterval of "
                    << fftWriteInterval_ << endl;
            }
            else
            {
                forAll(fOut, i)
                {
                    label freqI = (i + f0)*fftWriteInterval_;
                    fOut[i] = freq1[freqI];


                    PrmsfAve[i] = writeSurfaceData
                    (
                        outSurfDir,
                        fNameBase,
                        "Prmsf",
                        freq1[freqI],
                        surfPrmsf[i + f0],
                        procFaceAddr,
                        writePrmsf_
                    );

                    PSDfAve[i] = writeSurfaceData
                    (
                        outSurfDir,
                        fNameBase,
                        "PSDf",
                        freq1[freqI],
                        surfPSDf[i + f0],
                        procFaceAddr,
                        writePSDf_
                    );
                    writeSurfaceData
                    (
                        outSurfDir,
                        fNameBase,
                        "PSD",
                        freq1[freqI],
                        PSD(surfPSDf[i + f0]),
                        procFaceAddr,
                        writePSD_
                    );
                    writeSurfaceData
                    (
                        outSurfDir,
                        fNameBase,
                        "SPL",
                        freq1[freqI],
                        SPL(surfPSDf[i + f0]*deltaf, freq1[freqI]),
                        procFaceAddr,
                        writeSPL_
                    );
                }
            }

            if (Pstream::master())
            {
                {
                    auto filePtr = newFile(outDir/"Average_Prms_f");
                    auto& os = filePtr();

                    Info<< "    Writing " << os.relativeName() << endl;

                    writeFileHeader(os, "f [Hz]", "P(f) [Pa]", commonInfo);
                    writeFreqDataToFile(os, fOut, PrmsfAve);
                }
                {
                    auto filePtr = newFile(outDir/"Average_PSD_f_f");
                    auto& os = filePtr();

                    Info<< "    Writing " << os.relativeName() << endl;

                    writeFileHeader
                    (
                        os,
                        "f [Hz]",
                        "PSD(f) [PaPa_Hz]",
                        commonInfo
                    );
                    writeFreqDataToFile(os, fOut, PSDfAve);
                }
                {
                    auto filePtr = newFile(outDir/"Average_PSD_dB_Hz_f");
                    auto& os = filePtr();

                    Info<< "    Writing " << os.relativeName() << endl;

                    writeFileHeader
                    (
                        os,
                        "f [Hz]",
                        "PSD(f) [dB_Hz]",
                        commonInfo
                    );
                    writeFreqDataToFile(os, fOut, PSD(PSDfAve));
                }
                {
                    auto filePtr = newFile(outDir/"Average_SPL_dB_f");
                    auto& os = filePtr();

                    Info<< "    Writing " << os.relativeName() << endl;

                    writeFileHeader
                    (
                        os,
                        "f [Hz]",
                        "SPL(f) [dB]",
                        commonInfo
                    );
                    writeFreqDataToFile(os, fOut, SPL(PSDfAve*deltaf, fOut));
                }
            }
        }


        Info<< "Writing one-third octave surface data" << endl;
        {
            fileName outDir(outDirBase/"oneThirdOctave");
            fileName outSurfDir(filePath(outDir));

            scalarField PSDfAve(surfPrms13f.size(), Zero);
            scalarField Prms13fAve(surfPrms13f.size(), Zero);

            forAll(surfPrms13f, i)
            {
                writeSurfaceData
                (
                    outSurfDir,
                    fNameBase,
                    "SPL13",
                    octave13FreqCentre[i],
                    SPL(surfPrms13f[i], octave13FreqCentre[i]),
                    procFaceAddr,
                    writeOctaves_
                );

                Prms13fAve[i] =
                    surfaceAverage(surfPrms13f[i], procFaceAddr);
            }

            if (Pstream::master())
            {
                auto filePtr = newFile(outDir/"Average_SPL13_dB_fm");
                auto& os = filePtr();

                Info<< "    Writing " << os.relativeName() << endl;

                writeFileHeader
                (
                    os,
                    "fm [Hz]",
                    "SPL(fm) [dB]",
                    commonInfo
                );
                writeFreqDataToFile
                (
                    os,
                    octave13FreqCentre,
                    SPL(Prms13fAve, octave13FreqCentre)
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //
