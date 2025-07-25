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

#include "noiseModel.H"
#include "functionObject.H"
#include "argList.H"
#include "fft.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noiseModel, 0);
    defineRunTimeSelectionTable(noiseModel, dictionary);
}


const Foam::Enum<Foam::noiseModel::weightingType>
Foam::noiseModel::weightingTypeNames_
({
    {weightingType::none, "dB"},
    {weightingType::dBA, "dBA"},
    {weightingType::dBB, "dBB"},
    {weightingType::dBC, "dBC"},
    {weightingType::dBD, "dBD"},
});


void Foam::noiseModel::setOctaveBands
(
    const scalarField& f,
    const scalar fLower,
    const scalar fUpper,
    const scalar octave,
    labelList& fBandIDs,
    scalarField& fCentre
)
{
    // Set the indices of to the lower frequency bands for the input frequency
    // range. Ensure that the centre frequency passes though 1000 Hz

    // Low frequency bound given by:
    //     fLow = f0*(2^(0.5*bandI/octave))

    // Initial (lowest centre frequency)
    scalar fTest = 15.625;

    const scalar fRatio = pow(2, 1.0/octave);
    const scalar fRatioL2C = pow(2, 0.5/octave);

    // IDs of band IDs
    labelHashSet bandIDs(f.size());

    // Centre frequencies
    DynamicList<scalar> fc;
    DynamicList<scalar> missedBins;

    // Convert to lower band limit
    fTest /= fRatioL2C;
    while (fTest < fLower)
    {
        fTest *= fRatio;
    }

    forAll(f, i)
    {
        if (f[i] >= fTest)
        {
            // Advance band if appropriate
            label stepi = 0;
            while (f[i] > fTest)
            {
                if (stepi) missedBins.append(fTest/fRatio*fRatioL2C);
                fTest *= fRatio;
                ++stepi;
            }
            fTest /= fRatio;

            if (bandIDs.insert(i))
            {
                // Also store (next) centre frequency
                fc.append(fTest*fRatioL2C);
            }
            fTest *= fRatio;

            if (fTest > fUpper)
            {
                break;
            }
        }
    }

    fBandIDs = bandIDs.sortedToc();

    if (missedBins.size())
    {
        label nMiss = missedBins.size();
        label nTotal = nMiss + fc.size() - 1;
        WarningInFunction
            << "Empty bands found: " << nMiss << " of " << nTotal
            << " with centre frequencies " << flatOutput(missedBins)
            << endl;
    }

    if (fc.size())
    {
        // Remove the last centre frequency (beyond upper frequency limit)
        fc.pop_back();

        fCentre.transfer(fc);
    }
}


namespace Foam
{
    tmp<scalarField> safeLog10(const scalarField& fld)
    {
        auto tresult = tmp<scalarField>::New(fld.size(), -GREAT);
        auto& result = tresult.ref();

        forAll(result, i)
        {
            if (fld[i] > 0)
            {
                result[i] = log10(fld[i]);
            }
        }

        return tresult;
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Get bool option (eg, read/write) and provide Info feedback
static void readWriteOption
(
    const dictionary& dict,
    const word& lookup,
    bool& option
)
{
    dict.readIfPresent(lookup, option);

    Info<< "        " << lookup << ": " << (option ? "yes" : "no") << endl;
}

} // End namespace Foam


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::scalar Foam::noiseModel::checkUniformTimeStep
(
    const scalarList& times
) const
{
    scalar deltaT = 1;

    if (times.size() > 1)
    {
        // Uniform time step
        deltaT = (times.back() - times.front())/scalar(times.size() - 1);

        bool nonUniform = false;
        for (label i = 1; i < times.size() && !nonUniform; ++i)
        {
            const scalar dT = times[i] - times[i-1];

            nonUniform = (mag(dT/deltaT - 1) > 1e-8);
        }

        if (nonUniform)
        {
            WarningInFunction
                << "Detected non-uniform time step:"
                << " resulting FFT may be incorrect"
                << " (or the deltaT is just very small): " << deltaT
                << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unable to create FFT with 0 or 1 values"
            << exit(FatalError);
    }

    return deltaT;
}


bool Foam::noiseModel::validateBounds(const scalarList& p) const
{
    forAll(p, i)
    {
        if ((p[i] < minPressure_) || (p[i] > maxPressure_))
        {
            WarningInFunction
                << "Pressure data at position " << i
                << " is outside of permitted bounds:" << nl
                << "    pressure: " << p[i] << nl
                << "    minimum pressure: " << minPressure_ << nl
                << "    maximum pressure: " << maxPressure_ << nl
                << endl;

            return false;
        }
    }

    return true;
}


Foam::fileName Foam::noiseModel::baseFileDir(const label dataseti) const
{
    return
    (
        outputPrefix_
      / type()
      / ("input" + Foam::name(dataseti))
    );
}


void Foam::noiseModel::writeFileHeader
(
    Ostream& os,
    const string& x,
    const string& y,
    const UList<Tuple2<string, token>>& headerValues
) const
{
    writeHeader(os, x + " vs " + y);
    writeHeaderValue(os, "Lower frequency", fLower_);
    writeHeaderValue(os, "Upper frequency", fUpper_);
    writeHeaderValue(os, "Window model", windowModelPtr_->type());
    writeHeaderValue(os, "Window number", windowModelPtr_->nWindow());
    writeHeaderValue(os, "Window samples", windowModelPtr_->nSamples());
    writeHeaderValue(os, "Window overlap %", windowModelPtr_->overlapPercent());
    writeHeaderValue(os, "dBRef", dBRef_);

    for (const auto& hv : headerValues)
    {
        writeHeaderValue(os, hv.first(), hv.second());
    }

    writeCommented(os, x.substr(0, x.find(' ')));
    writeTabbed(os, y.substr(0, y.find(' ')));
    os  << nl;
}


void Foam::noiseModel::writeFreqDataToFile
(
    Ostream& os,
    const scalarField& f,
    const scalarField& fx
) const
{
    forAll(f, i)
    {
        if (f[i] >= fLower_ && f[i] <= fUpper_)
        {
            os  << f[i] << tab << fx[i] << nl;
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::uniformFrequencies
(
    const scalar deltaT,
    const bool check
) const
{
    const auto& window = windowModelPtr_();
    const label N = window.nSamples();

    auto tf = tmp<scalarField>::New(N/2 + 1, Zero);
    auto& f = tf.ref();

    const scalar deltaf = 1.0/(N*deltaT);

    label nFreq = 0;
    forAll(f, i)
    {
        f[i] = i*deltaf;

        if (f[i] > fLower_ && f[i] < fUpper_)
        {
            ++nFreq;
        }
    }

    if (check && nFreq == 0)
    {
        WarningInFunction
            << "No frequencies found in range "
            << fLower_ << " to " << fUpper_
            << endl;
    }

    return tf;
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::octaves
(
    const scalarField& data,
    const scalarField& f,
    const labelUList& freqBandIDs
) const
{
    if (freqBandIDs.size() < 2)
    {
        WarningInFunction
            << "Octave frequency bands are not defined "
            << "- skipping octaves calculation"
            << endl;

        return tmp<scalarField>::New();
    }

    auto toctData = tmp<scalarField>::New(freqBandIDs.size() - 1, Zero);
    auto& octData = toctData.ref();

    bitSet bandUsed(freqBandIDs.size() - 1);
    for (label bandI = 0; bandI < freqBandIDs.size() - 1; ++bandI)
    {
        label fb0 = freqBandIDs[bandI];
        label fb1 = freqBandIDs[bandI+1];

        if (fb0 == fb1) continue;

        for (label freqI = fb0; freqI < fb1; ++freqI)
        {
            label f0 = f[freqI];
            label f1 = f[freqI + 1];
            scalar dataAve = 0.5*(data[freqI] + data[freqI + 1]);
            octData[bandI] += dataAve*(f1 - f0);

            bandUsed.set(bandI);
        }
    }

    bandUsed.flip();
    labelList bandUnused = bandUsed.sortedToc();
    if (bandUnused.size())
    {
        WarningInFunction
            << "Empty bands found: " << bandUnused.size() << " of "
            << bandUsed.size() << endl;
    }

    return toctData;
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::Pf(const scalarField& p) const
{
    if (planInfo_.active)
    {
        if (p.size() != planInfo_.windowSize)
        {
            FatalErrorInFunction
                << "Expected pressure data to have " << planInfo_.windowSize
                << " values, but received " << p.size() << " values"
                << abort(FatalError);
        }

        List<double>& in = planInfo_.in;
        const List<double>& out = planInfo_.out;
        forAll(in, i)
        {
            in[i] = p[i];
        }

        ::fftw_execute(planInfo_.plan);

        const label n = planInfo_.windowSize;
        const label nBy2 = n/2;
        auto tresult = tmp<scalarField>::New(nBy2 + 1);
        auto& result = tresult.ref();

        // 0 th value = DC component
        // nBy2 th value = real only if n is even
        result[0] = out[0];
        for (label i = 1; i <= nBy2; ++i)
        {
            const auto re = out[i];
            const auto im = out[n - i];
            result[i] = sqrt(re*re + im*im);
        }

        return tresult;
    }

    return mag(fft::realTransform1D(p));
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::meanPf(const scalarField& p) const
{
    const auto& window = windowModelPtr_();
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    auto tmeanPf = tmp<scalarField>::New(N/2 + 1, Zero);
    auto& meanPf = tmeanPf.ref();

    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        meanPf += Pf(window.apply<scalar>(p, windowI));
    }

    meanPf /= scalar(nWindow);

    return tmeanPf;
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::RMSmeanPf
(
    const scalarField& p
) const
{
    const auto& window = windowModelPtr_();
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    scalarField RMSMeanPf(N/2 + 1, Zero);
    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        RMSMeanPf += sqr(Pf(window.apply<scalar>(p, windowI)));
    }

    return sqrt(RMSMeanPf/scalar(nWindow))/scalar(N);
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::PSDf
(
    const scalarField& p,
    const scalar deltaT
) const
{
    const auto& window = windowModelPtr_();
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    auto tpsd = tmp<scalarField>::New(N/2 + 1, Zero);
    auto& psd = tpsd.ref();

    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        psd += sqr(Pf(window.apply<scalar>(p, windowI)));
    }

    scalar fs =  1.0/deltaT;
    psd /= scalar(nWindow)*fs*N;

    // Scaling due to use of 1-sided FFT
    // Note: do not scale DC component
    psd *= 2;
    psd.first() /= 2;
    psd.last() /= 2;

    if (debug)
    {
        Pout<< "Average PSD: " << average(psd) << endl;
    }

    return tpsd;
}


Foam::scalar Foam::noiseModel::RAf(const scalar f) const
{
    const scalar c1 = sqr(12194.0);
    const scalar c2 = sqr(20.6);
    const scalar c3 = sqr(107.7);
    const scalar c4 = sqr(739.9);

    const scalar f2 = f*f;

    return
        c1*sqr(f2)
       /(
            (f2 + c2)*sqrt((f2 + c3)*(f2 + c4))*(f2 + c1)
        );
}


Foam::scalar Foam::noiseModel::gainA(const scalar f) const
{
    if (f < SMALL)
    {
        return 0;
    }

    return 20*log10(RAf(f)) - 20*log10(RAf(1000));
}


Foam::scalar Foam::noiseModel::RBf(const scalar f) const
{
    const scalar c1 = sqr(12194.0);
    const scalar c2 = sqr(20.6);
    const scalar c3 = sqr(158.5);

    const scalar f2 = f*f;

    return
        c1*f2*f
       /(
            (f2 + c2)*sqrt(f2 + c3)*(f2 + c1)
        );
}


Foam::scalar Foam::noiseModel::gainB(const scalar f) const
{
    if (f < SMALL)
    {
        return 0;
    }

    return 20*log10(RBf(f)) - 20*log10(RBf(1000));
}


Foam::scalar Foam::noiseModel::RCf(const scalar f) const
{
    const scalar c1 = sqr(12194.0);
    const scalar c2 = sqr(20.6);

    const scalar f2 = f*f;

    return c1*f2/((f2 + c2)*(f2 + c1));
}


Foam::scalar Foam::noiseModel::gainC(const scalar f) const
{
    if (f < SMALL)
    {
        return 0;
    }

    return 20*log10(RCf(f)) - 20*log10(RCf(1000));
}


Foam::scalar Foam::noiseModel::RDf(const scalar f) const
{
    const scalar f2 = f*f;

    const scalar hf =
        (sqr(1037918.48 - f2) + 1080768.16*f2)
       /(sqr(9837328 - f2) + 11723776*f2);

    return
        f/6.8966888496476e-5*Foam::sqrt(hf/((f2 + 79919.29)*(f2 + 1345600)));
}


Foam::scalar Foam::noiseModel::gainD(const scalar f) const
{
    if (f < SMALL)
    {
        return 0;
    }

    return 20*log10(RDf(f));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noiseModel::noiseModel
(
    const dictionary& dict,
    const objectRegistry& obr,
    const word& name,
    const bool readFields
)
:
    functionObjects::writeFile(obr, "noise"),
    dict_(dict),
    rhoRef_(1),
    nSamples_(65536),
    fLower_(25),
    fUpper_(10000),
    sampleFreq_(0),
    startTime_(0),
    windowModelPtr_(),
    SPLweighting_(weightingType::none),
    dBRef_(2e-5),
    minPressure_(-0.5*VGREAT),
    maxPressure_(0.5*VGREAT),
    outputPrefix_(),
    writePrmsf_(true),
    writeSPL_(true),
    writePSD_(true),
    writePSDf_(true),
    writeOctaves_(true),
    planInfo_()
{
    planInfo_.active = false;

    if (readFields)
    {
        read(dict);
    }

    if (debug)
    {
        writeWeightings();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::noiseModel::read(const dictionary& dict)
{
    if (!functionObjects::writeFile::read(dict))
    {
        return false;
    }

    // Re-assign defaults (to avoid stickiness)
    fLower_ = 25;
    fUpper_ = 10000;
    sampleFreq_ = 0;

    dict.readIfPresent("rhoRef", rhoRef_);
    dict.readIfPresent("N", nSamples_);
    dict.readIfPresentCompat("minFreq", {{"fl", 2312}}, fLower_);
    dict.readIfPresentCompat("maxFreq", {{"fu", 2312}}, fUpper_);
    dict.readIfPresent("sampleFreq", sampleFreq_);
    dict.readIfPresent("startTime", startTime_);
    dict.readIfPresent("minPressure", minPressure_);
    dict.readIfPresent("maxPressure", maxPressure_);
    dict.readIfPresent("outputPrefix", outputPrefix_);

    if (fLower_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Lower frequency bound must be greater than zero"
            << exit(FatalIOError);
    }

    if (fUpper_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Upper frequency bound must be greater than zero"
            << exit(FatalIOError);
    }

    if (fUpper_ < fLower_)
    {
        FatalIOErrorInFunction(dict)
            << "Upper frequency bound (" << fUpper_
            << ") must be greater than lower frequency bound ("
            << fLower_ << ")"
            << exit(FatalIOError);
    }

    Info<< "    Frequency bounds:" << nl
        << "        lower: " << fLower_ << nl
        << "        upper: " << fUpper_ << nl
        << "       sample: ";

    if (sampleFreq_ > 0)
    {
        Info<< sampleFreq_ << nl;
    }
    else
    {
        Info<< "auto" << nl;
    }

    weightingTypeNames_.readIfPresent("SPLweighting", dict, SPLweighting_);

    Info<< "    Weighting: " << weightingTypeNames_[SPLweighting_] << endl;

    if (dict.readIfPresent("dBRef", dBRef_))
    {
        Info<< "    Reference for dB calculation: " << dBRef_ << endl;
    }

    Info<< "    Write options:" << endl;
    dictionary optDict(dict.subOrEmptyDict("writeOptions"));
    readWriteOption(optDict, "writePrmsf", writePrmsf_);
    readWriteOption(optDict, "writeSPL", writeSPL_);
    readWriteOption(optDict, "writePSD", writePSD_);
    readWriteOption(optDict, "writePSDf", writePSDf_);
    readWriteOption(optDict, "writeOctaves", writeOctaves_);


    windowModelPtr_ = windowModel::New(dict, nSamples_);

    cleanFFTW();

    const label windowSize = windowModelPtr_->nSamples();

    if (windowSize > 1)
    {
        planInfo_.active = true;
        planInfo_.windowSize = windowSize;
        planInfo_.in.setSize(windowSize);
        planInfo_.out.setSize(windowSize);

        // Using real to half-complex fftw 'kind'
        planInfo_.plan =
            fftw_plan_r2r_1d
            (
                windowSize,
                planInfo_.in.data(),
                planInfo_.out.data(),
                FFTW_R2HC,
                windowSize <= 8192 ? FFTW_MEASURE : FFTW_ESTIMATE
            );
    }

    Info<< endl;

    return true;
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::PSD
(
    const scalarField& PSDf
) const
{
    return 10*safeLog10(PSDf/sqr(dBRef_));
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::SPL
(
    const scalarField& Prms2,
    const scalar f
) const
{
    tmp<scalarField> tspl(10*safeLog10(Prms2/sqr(dBRef_)));
    scalarField& spl = tspl.ref();

    switch (SPLweighting_)
    {
        case weightingType::none:
        {
            break;
        }
        case weightingType::dBA:
        {
            spl += gainA(f);
            break;
        }
        case weightingType::dBB:
        {
            spl += gainB(f);
            break;
        }
        case weightingType::dBC:
        {
            spl += gainC(f);
            break;
        }
        case weightingType::dBD:
        {
            spl += gainD(f);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown weighting " << weightingTypeNames_[SPLweighting_]
                << abort(FatalError);
        }
    }

    return tspl;
}


Foam::tmp<Foam::scalarField> Foam::noiseModel::SPL
(
    const scalarField& Prms2,
    const scalarField& f
) const
{
    tmp<scalarField> tspl(10*safeLog10(Prms2/sqr(dBRef_)));
    scalarField& spl = tspl.ref();

    switch (SPLweighting_)
    {
        case weightingType::none:
        {
            break;
        }
        case weightingType::dBA:
        {
            forAll(spl, i)
            {
                spl[i] += gainA(f[i]);
            }
            break;
        }
        case weightingType::dBB:
        {
            forAll(spl, i)
            {
                spl[i] += gainB(f[i]);
            }
            break;
        }
        case weightingType::dBC:
        {
            forAll(spl, i)
            {
                spl[i] += gainC(f[i]);
            }
            break;
        }
        case weightingType::dBD:
        {
            forAll(spl, i)
            {
                spl[i] += gainD(f[i]);
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown weighting " << weightingTypeNames_[SPLweighting_]
                << abort(FatalError);
        }
    }

    return tspl;
}


void Foam::noiseModel::cleanFFTW()
{
    if (planInfo_.active)
    {
        planInfo_.active = false;
        fftw_destroy_plan(planInfo_.plan);
        fftw_cleanup();
    }
}


void Foam::noiseModel::writeWeightings() const
{
    scalar f0 = 10;
    scalar f1 = 20000;

    OFstream osA("noiseModel-weight-A");
    OFstream osB("noiseModel-weight-B");
    OFstream osC("noiseModel-weight-C");
    OFstream osD("noiseModel-weight-D");

    for (label f = f0; f <= f1; ++f)
    {
        osA << f << " " << gainA(f) << nl;
        osB << f << " " << gainB(f) << nl;
        osC << f << " " << gainC(f) << nl;
        osD << f << " " << gainD(f) << nl;
    }
}


// ************************************************************************* //
