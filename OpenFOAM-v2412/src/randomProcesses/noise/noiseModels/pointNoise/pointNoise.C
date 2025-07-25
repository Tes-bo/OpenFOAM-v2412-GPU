/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "pointNoise.H"
#include "argList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace noiseModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointNoise, 0);
addToRunTimeSelectionTable(noiseModel, pointNoise, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void pointNoise::filterTimeData
(
    const scalarField& t0,
    const scalarField& p0,
    scalarField& t,
    scalarField& p
) const
{
    DynamicList<scalar> tf(t0.size());
    DynamicList<scalar> pf(t0.size());

    forAll(t0, timeI)
    {
        if (t0[timeI] >= startTime_)
        {
            tf.append(t0[timeI]);
            pf.append(p0[timeI]);
        }
    }

    t.transfer(tf);
    p.transfer(pf);
}


void pointNoise::processData
(
    const label dataseti,
    const Function1Types::CSV<scalar>& data
)
{
    if (!Pstream::master())
    {
        // Only ever called on master, report if we have odd logic...
        WarningInFunction
            << "Currently only to be called from master process..." << endl;
        return;
    }

    Info<< "Reading data file: "
        << fileObr_.time().relativePath(data.fName()) << endl;

    const word fNameBase(data.fName().stem());

    // Time and pressure history data
    scalarField t, p;
    filterTimeData(data.x(), data.y(), t, p);

    Info<< "    read " << t.size() << " values" << nl << endl;

    if (!validateBounds(p))
    {
        Info<< "No noise data generated" << endl;
        return;
    }

    Info<< "Creating noise FFT" << endl;

    const scalar deltaT =
    (
        sampleFreq_ > 0
      ? (1.0/sampleFreq_)
      : checkUniformTimeStep(t)
    );

    // Apply conversions
    p *= rhoRef_;
    p -= average(p);

    // Determine the windowing
    windowModelPtr_->validate(t.size());
    const windowModel& win = windowModelPtr_();
    const scalar deltaf = 1.0/(deltaT*win.nSamples());
    const fileName outDir(baseFileDir(dataseti)/fNameBase);


    // Narrow band data
    // ----------------

    scalarField f(uniformFrequencies(deltaT, true));

    // RMS pressure [Pa]
    if (writePrmsf_)
    {
        auto filePtr = newFile(outDir/"Prms_f");
        auto& os = filePtr();

        Info<< "    Writing " << os.name() << endl;

        writeFileHeader(os, "f [Hz]", "Prms(f) [Pa]");
        writeFreqDataToFile(os, f, RMSmeanPf(p));
    }

    // PSD [Pa^2/Hz]
    const scalarField PSDf(this->PSDf(p, deltaT));

    if (writePSDf_)
    {
        auto filePtr = newFile(outDir/"PSD_f");
        auto& os = filePtr();

        Info<< "    Writing " << os.relativeName() << endl;

        writeFileHeader(os, "f [Hz]", "PSD(f) [PaPa_Hz]");
        writeFreqDataToFile(os, f, PSDf);
    }

    // PSD [dB/Hz]
    if (writePSD_)
    {
        auto filePtr = newFile(outDir/"PSD_dB_Hz_f");
        auto& os = filePtr();

        Info<< "    Writing " << os.relativeName() << endl;

        writeFileHeader(os, "f [Hz]", "PSD(f) [dB_Hz]");
        writeFreqDataToFile(os, f, PSD(PSDf));
    }

    // SPL [dB]
    if (writeSPL_)
    {
        auto filePtr = newFile(outDir/"SPL_dB_f");
        auto& os = filePtr();

        Info<< "    Writing " << os.relativeName() << endl;

        writeFileHeader
        (
            os,
            "f [Hz]",
            "SPL(f) [" + weightingTypeNames_[SPLweighting_] + "]"
        );
        writeFreqDataToFile(os, f, SPL(PSDf*deltaf, f));
    }

    if (writeOctaves_)
    {
        labelList octave13BandIDs;
        scalarField octave13FreqCentre;
        setOctaveBands
        (
            f,
            fLower_,
            fUpper_,
            3,
            octave13BandIDs,
            octave13FreqCentre
        );


        // 1/3 octave data
        // ---------------

        // Integrated PSD = P(rms)^2 [Pa^2]
        scalarField Prms13f(octaves(PSDf, f, octave13BandIDs));

        auto filePtr = newFile(outDir/"SPL13_dB_fm");
        auto& os = filePtr();

        Info<< "    Writing " << os.relativeName() << endl;

        writeFileHeader
        (
            os,
            "fm [Hz]",
            "SPL(fm) [" + weightingTypeNames_[SPLweighting_] + "]"
        );
        writeFreqDataToFile
        (
            os,
            octave13FreqCentre,
            SPL(Prms13f, octave13FreqCentre)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointNoise::pointNoise
(
    const dictionary& dict,
    const objectRegistry& obr,
    const word& name,
    const bool readFields
)
:
    noiseModel(dict, obr, name, false)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointNoise::calculate()
{
    // Point data only handled by master
    if (!Pstream::master())
    {
        return;
    }

    forAll(inputFileNames_, filei)
    {
        fileName fName = inputFileNames_[filei];
        fName.expand();

        if (!fName.isAbsolute())
        {
            fName = argList::envGlobalPath()/fName;
        }

        Function1Types::CSV<scalar> data("pressure", dict_, nullptr, fName);
        processData(filei, data);
    }
}


bool pointNoise::read(const dictionary& dict)
{
    if (noiseModel::read(dict))
    {
        if (!dict.readIfPresent("files", inputFileNames_))
        {
            inputFileNames_.resize(1);

            // Note: lookup uses same keyword as used by the CSV constructor
            dict.readEntry("file", inputFileNames_.first());
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //
