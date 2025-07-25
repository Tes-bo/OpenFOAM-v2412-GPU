/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
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

#include "streamFunctionWaveModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(streamFunction, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        streamFunction,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::streamFunction::eta
(
    const scalar h,
    const scalar kx,
    const scalar ky,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase
) const
{

   const scalar k = sqrt(kx*kx + ky*ky);
   scalar strfnAux = 0.0;
   forAll(Ejs_, iterSF)
   {
        strfnAux +=
            Ejs_[iterSF]*cos((iterSF + 1) *(kx*x + ky*y - omega*t + phase));
   }

   return (1/k)*strfnAux;
}

Foam::vector Foam::waveModels::streamFunction::Uf
(
    const scalar h,
    const scalar kx,
    const scalar ky,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase,
    const scalar z
) const
{
    const scalar k = sqrt(kx*kx + ky*ky);
    const scalar phaseTot = kx*x + ky*y - omega*t + phase;

    scalar u = 0.0;
    scalar w = 0.0;

    forAll(Bjs_, iterSF2)
    {
        u +=
            (iterSF2 + 1)*Bjs_[iterSF2]*cosh((iterSF2 + 1)*k*z)
           /cosh((iterSF2 + 1)*k*h)*cos((iterSF2 + 1)*phaseTot);

        w +=
            (iterSF2 + 1)*Bjs_[iterSF2]*sinh((iterSF2 + 1)*k*z)
           /cosh((iterSF2 + 1)*k*h)*sin((iterSF2 + 1)*phaseTot);
    }

    u = waveLength_/T - uMean_ + sqrt(mag(g_)/k)*u;
    w = sqrt(mag(g_)/k)*w;

    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::waveModels::streamFunction::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    const scalar waveOmega = mathematical::twoPi/wavePeriod_;
    const scalar waveK = mathematical::twoPi/waveLength_;
    const scalar waveKx = waveK*cos(waveAngle_);
    const scalar waveKy = waveK*sin(waveAngle_);

    forAll(level, paddlei)
    {
        const scalar eta =
            this->eta
            (
                waterDepthRef_,
                waveKx,
                waveKy,
                wavePeriod_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveOmega,
                t,
                wavePhase_
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


void Foam::waveModels::streamFunction::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
    const scalar waveOmega = mathematical::twoPi/wavePeriod_;
    const scalar waveK = mathematical::twoPi/waveLength_;
    const scalar waveKx = waveK*cos(waveAngle_);
    const scalar waveKy = waveK*sin(waveAngle_);

    forAll(U_, facei)
    {
        // Fraction of geometry represented by paddle - to be set
        scalar fraction = 1;

        // Height - to be set
        scalar z = 0;

        setPaddlePropeties(level, facei, fraction, z);

        if (fraction > 0)
        {
            const label paddlei = faceToPaddle_[facei];

            const vector Uf = this->Uf
            (
                waterDepthRef_,
                waveKx,
                waveKy,
                wavePeriod_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveOmega,
                t,
                wavePhase_,
                z
            );

            U_[facei] = fraction*Uf*tCoeff;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::streamFunction::streamFunction
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    regularWaveModel(dict, mesh, patch, false),
    uMean_(0),
    Bjs_(),
    Ejs_()
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::streamFunction::readDict(const dictionary& overrideDict)
{
    if (regularWaveModel::readDict(overrideDict))
    {
        overrideDict.readEntry("uMean", uMean_);
        overrideDict.readEntry("waveLength", waveLength_);
        overrideDict.readEntry("Bjs", Bjs_);
        overrideDict.readEntry("Ejs", Ejs_);

        return true;
    }

    return false;
}


void Foam::waveModels::streamFunction::info(Ostream& os) const
{
    regularWaveModel::info(os);

    os  << "    uMean : " << uMean_ << nl
        << "    Stream function wavelength : " << waveLength_ << nl
        << "    Bj coefficients : " << Bjs_ << nl
        << "    Ej coefficients : " << Ejs_ << nl;
}


// ************************************************************************* //
