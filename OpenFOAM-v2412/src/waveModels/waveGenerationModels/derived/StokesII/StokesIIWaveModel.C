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

#include "StokesIIWaveModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(StokesII, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        StokesII,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::StokesII::eta
(
    const scalar H,
    const scalar h,
    const scalar Kx,
    const scalar x,
    const scalar Ky,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase
) const
{
    const scalar k = sqrt(Kx*Kx + Ky*Ky);
    const scalar sigma = tanh(k*h);
    const scalar phaseTot = Kx*x + Ky*y - omega*t + phase;

    return
        H*0.5*cos(phaseTot)
      + k*H*H/4.0*(3.0 - sigma*sigma)/(4.0*pow3(sigma))*cos(2.0*phaseTot);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::vector Foam::waveModels::StokesII::UfBase
(
    const scalar H,
    const scalar h,
    const scalar Kx,
    const scalar x,
    const scalar Ky,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase,
    const scalar z
) const
{
    const scalar k = sqrt(Kx*Kx + Ky*Ky);
    const scalar phaseTot = Kx*x + Ky*y - omega*t + phase;

    scalar u =
        H*0.5*omega*cos(phaseTot)*cosh(k*z)/sinh(k*h)
      + 3.0/4.0*H*H/4.0*omega*k*cosh(2.0*k*z)/pow4(sinh(k*h))*cos(2.0*phaseTot);
    scalar w =
        H*0.5*omega*sin(phaseTot)*sinh(k*z)/sinh(k*h)
      + 3.0/4.0*H*H/4.0*omega*k*sinh(2.0*k*z)/pow4(sinh(k*h))*sin(2.0*phaseTot);
    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


void Foam::waveModels::StokesII::setLevel
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
                waveHeight_,
                waterDepthRef_,
                waveKx,
                xPaddle_[paddlei],
                waveKy,
                yPaddle_[paddlei],
                waveOmega,
                t,
                wavePhase_
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::StokesII::StokesII
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    StokesI(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::StokesII::readDict(const dictionary& overrideDict)
{
    if (StokesI::readDict(overrideDict))
    {
        return true;
    }

    return false;
}


void Foam::waveModels::StokesII::info(Ostream& os) const
{
    StokesI::info(os);
}


// ************************************************************************* //
