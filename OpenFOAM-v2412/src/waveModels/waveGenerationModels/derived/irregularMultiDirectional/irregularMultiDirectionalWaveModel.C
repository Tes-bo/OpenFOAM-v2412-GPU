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

#include "irregularMultiDirectionalWaveModel.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(irregularMultiDirectional, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        irregularMultiDirectional,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::irregularMultiDirectional::eta
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

Foam::scalar Foam::waveModels::irregularMultiDirectional::waveLength
(
    const scalar h,
    const scalar T
) const
{
    scalar L0 = mag(g_)*T*T/(2.0*mathematical::pi);
    scalar L = L0;

    for (label i=1; i<=100; ++i)
    {
        L = L0*tanh(2.0*mathematical::pi*h/L);
    }

    return L;
}


Foam::vector Foam::waveModels::irregularMultiDirectional::Uf
(
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar t,
    const scalar z
) const
{
    scalar u = 0.0;
    scalar v = 0.0;
    scalar w = 0.0;

    forAll(irregWaveHeights_, ii)
    {
        forAll(irregWaveHeights_[ii], jj)
        {
            scalar waveKs = mathematical::twoPi/irregWaveLengths_[ii][jj];
            scalar waveOmegas = mathematical::twoPi/irregWavePeriods_[ii][jj];

            scalar phaseTot =
                waveKs*x*cos(irregWaveDirs_[ii][jj])
              + waveKs*y*sin(irregWaveDirs_[ii][jj])
              - waveOmegas*t
              + irregWavePhases_[ii][jj];

            const vector Uf = uMultiDirec
            (
                irregWaveHeights_[ii][jj],
                waveOmegas,
                phaseTot,
                waveKs,
                z,
                h,
                irregWaveDirs_[ii][jj]
            );

            u += Uf[0];
            v += Uf[1];
            w += Uf[2];
        }
    }

    return vector(u, v, w);
}


void Foam::waveModels::irregularMultiDirectional::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    forAll(level, paddlei)
    {
        scalar eta = 0;

        forAll(irregWaveHeights_, ii)
        {
            forAll(irregWaveHeights_[ii], jj)
            {
                scalar waveKs = mathematical::twoPi/irregWaveLengths_[ii][jj];
                scalar waveOmegas =
                    mathematical::twoPi/irregWavePeriods_[ii][jj];

                eta +=
                    this->eta
                    (
                        irregWaveHeights_[ii][jj],
                        waveKs*cos(irregWaveDirs_[ii][jj]),
                        xPaddle_[paddlei],
                        waveKs*sin(irregWaveDirs_[ii][jj]),
                        yPaddle_[paddlei],
                        waveOmegas,
                        t,
                        irregWavePhases_[ii][jj]
                    );
            }
        }

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


Foam::vector Foam::waveModels::irregularMultiDirectional::uMultiDirec
(
    const scalar irregH,
    const scalar irregWaveOmega,
    const scalar pha,
    const scalar irregWaveKs,
    const scalar zz,
    const scalar hh,
    const scalar irregDir
) const
{
    const scalar ksh = irregWaveKs*hh;
    const scalar ksz = irregWaveKs*zz;

    scalar u =
        irregH*0.5*irregWaveOmega*cos(pha)*(cosh(ksz)/sinh(ksh))*cos(irregDir);

    scalar v =
        irregH*0.5*irregWaveOmega*cos(pha)*(cosh(ksz)/sinh(ksh))*sin(irregDir);

    scalar w =
        irregH*0.5*irregWaveOmega*sin(pha)*(sinh(ksz)/sinh(ksh));

    return vector(u, v, w);
}


void Foam::waveModels::irregularMultiDirectional::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
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
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                t,
                z
            );

            U_[facei] = fraction*Uf*tCoeff;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::irregularMultiDirectional::irregularMultiDirectional
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    irregularWaveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::irregularMultiDirectional::readDict
(
    const dictionary& overrideDict
)
{
    if (irregularWaveModel::readDict(overrideDict))
    {
        readEntry("wavePeriods", irregWavePeriods_);
        readEntry("waveHeights", irregWaveHeights_);
        readEntry("wavePhases", irregWavePhases_);
        readEntry("waveDirs", irregWaveDirs_);

        irregWaveLengths_ = irregWaveHeights_;

        forAll(irregWaveHeights_, ii)
        {
            forAll(irregWaveHeights_[ii], jj)
            {
                irregWaveLengths_[ii][jj] =
                    waveLength(waterDepthRef_, irregWavePeriods_[ii][jj]);
                irregWaveDirs_[ii][jj] =
                    degToRad(irregWaveDirs_[ii][jj]);
            }
        }

        return true;
    }

    return false;
}


void Foam::waveModels::irregularMultiDirectional::info(Ostream& os) const
{
    irregularWaveModel::info(os);

    os  << "    Wave periods    : " << irregWavePeriods_.size() << nl
        << "    Wave heights    : " << irregWaveHeights_.size() << nl
        << "    Wave phases     : " << irregWavePhases_.size() << nl
        << "    Wave lengths    : " << irregWaveLengths_.size() << nl
        << "    Wave directions : " << irregWaveDirs_.size() << nl;
}


// ************************************************************************* //
