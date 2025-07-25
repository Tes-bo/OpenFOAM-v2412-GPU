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

#include "StokesIWaveModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(StokesI, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        StokesI,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::StokesI::eta
(
    const scalar H,
    const scalar Kx,
    const scalar x,
    const scalar Ky,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase
) const
{
    scalar phaseTot = Kx*x + Ky*y - omega*t + phase;

    return H*0.5*cos(phaseTot);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::StokesI::waveLength
(
    const scalar h,
    const scalar T
) const
{
    scalar L0 = mag(g_)*T*T/(2.0*mathematical::pi);
    scalar L = L0;

    for (int i=1; i<=100; i++)
    {
        L = L0*tanh(2.0*mathematical::pi*h/L);
    }

    return L;
}


Foam::vector Foam::waveModels::StokesI::UfBase
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
    scalar k = sqrt(Kx*Kx + Ky*Ky);
    scalar phaseTot = Kx*x + Ky*y - omega*t + phase;

    scalar u = H*0.5*omega*cos(phaseTot)*cosh(k*z)/sinh(k*h);
    scalar w = H*0.5*omega*sin(phaseTot)*sinh(k*z)/sinh(k*h);
    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


void Foam::waveModels::StokesI::setLevel
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


void Foam::waveModels::StokesI::setVelocity
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

            const vector Uf = UfBase
            (
                waveHeight_,
                waterDepthRef_,
                waveKx,
                xPaddle_[paddlei],
                waveKy,
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

Foam::waveModels::StokesI::StokesI
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    regularWaveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::StokesI::readDict(const dictionary& overrideDict)
{
    if (regularWaveModel::readDict(overrideDict))
    {
        waveLength_ = waveLength(waterDepthRef_, wavePeriod_);

        return true;
    }

    return false;
}


void Foam::waveModels::StokesI::info(Ostream& os) const
{
    regularWaveModel::info(os);

    os  << "    Wave type: " << waveType() << nl;
}


// ************************************************************************* //
