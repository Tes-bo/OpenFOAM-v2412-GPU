/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "noiseFFT.H"
#include "IFstream.H"
#include "DynamicList.H"
#include "mathematicalConstants.H"
#include "HashSet.H"
#include "fft.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * //

Foam::scalar Foam::noiseFFT::p0 = 2e-5;


Foam::tmp<Foam::scalarField> Foam::noiseFFT::frequencies
(
    const label N,
    const scalar deltaT
)
{
    auto tf = tmp<scalarField>::New(N/2, Zero);
    auto& f = tf.ref();

    const scalar deltaf = 1.0/(N*deltaT);
    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return tf;
}


Foam::tmp<Foam::scalarField> Foam::noiseFFT::PSD(const scalarField& PSDf)
{
    return 10*log10(PSDf/sqr(p0));
}


Foam::tmp<Foam::scalarField> Foam::noiseFFT::SPL(const scalarField& Prms2)
{
    return 10*log10(Prms2/sqr(p0));
}


void Foam::noiseFFT::octaveBandInfo
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

    // Convert to lower band limit
    fTest /= fRatioL2C;

    forAll(f, i)
    {
        if (f[i] >= fTest)
        {
            // Advance band if appropriate
            while (f[i] > fTest)
            {
                fTest *= fRatio;
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

    if (fc.size())
    {
        // Remove the last centre frequency (beyond upper frequency limit)
        fc.pop_back();

        fCentre.transfer(fc);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noiseFFT::noiseFFT(const scalar deltaT, const label windowSize)
:
    scalarField(),
    deltaT_(deltaT)
{
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
    else
    {
        planInfo_.active = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::noiseFFT::~noiseFFT()
{
    if (planInfo_.active)
    {
        planInfo_.active = false;
        fftw_destroy_plan(planInfo_.plan);
        fftw_cleanup();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noiseFFT::setData(scalarList& data)
{
    this->transfer(data);

    scalarField& p = *this;
    p -= average(p);
}


void Foam::noiseFFT::setData(const fileName& pFileName, const label skip)
{
    // Construct pressure data file
    IFstream pFile(pFileName);

    // Check pFile stream is OK
    if (!pFile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << pFileName
            << exit(FatalError);
    }

    if (skip)
    {
        scalar dummyt, dummyp;

        for (label i = 0; i < skip; ++i)
        {
            pFile >> dummyt;

            if (!pFile.good() || pFile.eof())
            {
                FatalErrorInFunction
                    << "Number of points in file " << pFileName
                    << " is less than the number to be skipped = " << skip
                    << exit(FatalError);
            }

            pFile >> dummyp;
        }
    }

    scalar t = 0, T0 = 0, T1 = 0;
    DynamicList<scalar> pData(1024);
    label i = 0;

    while (!(pFile >> t).eof())
    {
        if (i == 0)
        {
            T0 = t;
        }

        T1 = t;
        pFile >> pData(i++);
    }

    // Note: assumes fixed time step
    deltaT_ = (T1 - T0)/pData.size();

    this->transfer(pData);

    scalarField& p = *this;
    p -= average(p);
}


Foam::graph Foam::noiseFFT::pt() const
{
    scalarField t(size());
    forAll(t, i)
    {
        t[i] = i*deltaT_;
    }

    return graph
    (
        "p(t)",
        "t [s]",
        "p(t) [Pa]",
        t,
        *this
    );
}


Foam::tmp<Foam::scalarField> Foam::noiseFFT::Pf
(
    const tmp<scalarField>& tpn
) const
{
    if (planInfo_.active)
    {
        const scalarField& pn = tpn();
        if (pn.size() != planInfo_.windowSize)
        {
            FatalErrorInFunction
                << "Expected pressure data to have " << planInfo_.windowSize
                << " values, but received " << pn.size() << " values"
                << abort(FatalError);
        }

        List<double>& in = planInfo_.in;
        const List<double>& out = planInfo_.out;
        forAll(in, i)
        {
            in[i] = pn[i];
        }
        tpn.clear();

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

    return mag(fft::realTransform1D(tpn));
}


Foam::graph Foam::noiseFFT::meanPf(const windowModel& window) const
{
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    scalarField meanPf(N/2 + 1, Zero);

    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        meanPf += Pf(window.apply<scalar>(*this, windowI));
    }

    meanPf /= scalar(nWindow);

    scalar deltaf = 1.0/(N*deltaT_);
    scalarField f(meanPf.size());
    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return graph
    (
        "P(f)",
        "f [Hz]",
        "P(f) [Pa]",
        f,
        meanPf
    );
}


Foam::graph Foam::noiseFFT::RMSmeanPf(const windowModel& window) const
{
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    scalarField RMSMeanPf(N/2 + 1, Zero);
    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        RMSMeanPf += sqr(Pf(window.apply<scalar>(*this, windowI)));
    }

    RMSMeanPf = sqrt(RMSMeanPf/scalar(nWindow))/scalar(N);

    scalar deltaf = 1.0/(N*deltaT_);
    scalarField f(RMSMeanPf.size());
    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return graph
    (
        "Prms(f)",
        "f [Hz]",
        "Prms(f) [Pa]",
        f,
        RMSMeanPf
    );
}


Foam::graph Foam::noiseFFT::PSDf(const windowModel& window) const
{
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    scalarField psd(N/2 + 1, Zero);

    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        psd += sqr(Pf(window.apply<scalar>(*this, windowI)));
    }

    scalar deltaf = 1.0/(N*deltaT_);
    scalar fs =  1.0/deltaT_;
    psd /= scalar(nWindow)*fs*N;

    // Scaling due to use of 1-sided FFT
    // Note: do not scale DC component
    psd *= 2;
    psd.first() /= 2;
    psd.last() /= 2;

    scalarField f(psd.size());

    if (0) // if (debug)
    {
        Pout<< "Average PSD: " << average(psd) << endl;
    }

    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return graph
    (
        "PSD(f)",
        "f [Hz]",
        "PSD(f) [PaPa_Hz]",
        f,
        psd
    );
}


Foam::graph Foam::noiseFFT::PSD(const graph& gPSDf) const
{
    return graph
    (
        "PSD(f)",
        "f [Hz]",
        "PSD_dB(f) [dB_Hz]",
        gPSDf.x(),
        10*log10(gPSDf.y()/sqr(p0))
    );
}


Foam::graph Foam::noiseFFT::octaves
(
    const graph& g,
    const labelUList& freqBandIDs
) const
{
    if (freqBandIDs.size() < 2)
    {
        WarningInFunction
            << "Octave frequency bands are not defined "
            << "- skipping octaves calculation"
            << endl;

        return graph
        (
            "octave",
            "x",
            "y",
            scalarField(),
            scalarField()
        );
    }

    const scalarField& f = g.x();
    const scalarField& data = g.y();

    scalarField octData(freqBandIDs.size() - 1, Zero);
    scalarField fm(freqBandIDs.size() -1, Zero);

    for (label bandI = 0; bandI < freqBandIDs.size() - 1; ++bandI)
    {
        label fb0 = freqBandIDs[bandI];
        label fb1 = freqBandIDs[bandI+1];
        fm[bandI] = f[fb0];

        if (fb0 == fb1) continue;

        for (label freqI = fb0; freqI < fb1; ++freqI)
        {
            label f0 = f[freqI];
            label f1 = f[freqI + 1];
            scalar dataAve = 0.5*(data[freqI] + data[freqI + 1]);
            octData[bandI] += dataAve*(f1 - f0);
        }
    }

    return graph
    (
        "octaves(f)",
        "fm [Hz]",
        "octave data",
        fm,
        octData
    );
}


Foam::scalar Foam::noiseFFT::dbToPa(const scalar db) const
{
    return p0*pow(10.0, db/20.0);
}


Foam::tmp<Foam::scalarField> Foam::noiseFFT::dbToPa
(
    const tmp<scalarField>& db
) const
{
    return p0*pow(10.0, db/20.0);
}


// ************************************************************************* //
