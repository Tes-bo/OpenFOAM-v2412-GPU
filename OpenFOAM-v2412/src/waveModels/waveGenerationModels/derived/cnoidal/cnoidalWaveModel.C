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

#include "cnoidalWaveModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "Elliptic.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(cnoidal, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        cnoidal,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::waveModels::cnoidal::initialise
(
    const scalar H,
    const scalar d,
    const scalar T,
    scalar& mOut,
    scalar& LOut
) const
{
    const scalar mTolerance = 0.0001;
    scalar mElliptic = 0.5;
    scalar LElliptic = 0;
    scalar phaseSpeed = 0;

    scalar mError = 0.0;
    scalar mMinError = GREAT;

    while (mElliptic < 1.0)
    {
        scalar KElliptic, EElliptic;
        Elliptic::ellipticIntegralsKE(mElliptic, KElliptic, EElliptic);

        LElliptic = KElliptic*sqrt(16.0*pow3(d)*mElliptic/(3.0*H));

        phaseSpeed =
            sqrt(mag(g_)*d)
           *(1.0 - H/d/2.0 + H/d/mElliptic*(1.0 - 3.0/2.0*EElliptic/KElliptic));

        mError = mag(T - LElliptic/phaseSpeed);

        if (mError <= mMinError)
        {
            mOut = mElliptic;
            LOut = LElliptic;
            mMinError = mError;
        }

        mElliptic += mTolerance;
    }
}


Foam::scalar Foam::waveModels::cnoidal::eta
(
    const scalar H,
    const scalar m,
    const scalar kx,
    const scalar ky,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar t
) const
{
    scalar K, E;
    Elliptic::ellipticIntegralsKE(m, K, E);

    const scalar uCnoidal =
        K/mathematical::pi*(kx*x + ky*y - 2.0*mathematical::pi*t/T);

    scalar sn, cn, dn;
    Elliptic::JacobiSnCnDn(uCnoidal, m, sn, cn, dn);

    return H*((1.0 - E/K)/m - 1.0 + sqr(cn));
}


Foam::scalar Foam::waveModels::cnoidal::eta1D
(
    const scalar H,
    const scalar m,
    const scalar t,
    const scalar T
) const
{
    scalar K, E;
    Elliptic::ellipticIntegralsKE(m, K, E);

    const scalar uCnoidal = -2.0*K*(t/T);

    scalar sn, cn, dn;
    Elliptic::JacobiSnCnDn(uCnoidal, m, sn, cn, dn);

    return H*((1.0 - E/K)/m - 1.0 + sqr(cn));
}


Foam::scalar Foam::waveModels::cnoidal::etaMeanSq
(
    const scalar H,
    const scalar m,
    const scalar T
) const
{
    scalar eta = 0;
    scalar etaSumSq = 0;

    for (int i=0; i<1000; i++)
    {
        eta = eta1D(H, m, i*T/(1000.0), T);
        etaSumSq += eta*eta;
    }

    etaSumSq /= 1000.0;
    return etaSumSq;
}


Foam::vector Foam::waveModels::cnoidal::dEtaDx
(
    const scalar H,
    const scalar m,
    const scalar uCnoidal,
    const scalar L,
    const scalar K,
    const scalar E
) const
{
    const scalar dudx = 2.0*K/L;
    const scalar dudxx = 2.0*K/L*dudx;
    const scalar dudxxx = 2.0*K/L*dudxx;

    scalar sn, cn, dn;
    Elliptic::JacobiSnCnDn(uCnoidal, m, sn, cn, dn);

    scalar d1 = -2.0*H*cn*dn*sn*dudx;
    scalar d2 = 2.0*H*(dn*dn*sn*sn - cn*cn*dn*dn + m*cn*cn*sn*sn)*dudxx;
    scalar d3 =
        8.0*H
       *(
            cn*sn*dn*dn*dn*(-4.0 - 2.0*m)
          + 4.0*m*cn*sn*sn*sn*dn
          - 2.0*m*cn*cn*cn*sn*dn
        )
       *dudxxx;

    return vector(d1, d2, d3);
}


Foam::vector Foam::waveModels::cnoidal::Uf
(
    const scalar H,
    const scalar h,
    const scalar m,
    const scalar kx,
    const scalar ky,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar t,
    const scalar z
) const
{
    scalar K, E;
    Elliptic::ellipticIntegralsKE(m, K, E);

    const scalar uCnoidal =
        K/mathematical::pi*(kx*x + ky*y - 2.0*mathematical::pi*t/T);
    const scalar k = sqrt(kx*kx + ky*ky);
    const scalar L = 2.0*mathematical::pi/k;
    const scalar c = L/T;

    const scalar etaCN = eta(H, m, kx, ky, T, x, y, t);
    const vector etaX = this->dEtaDx(H, m, uCnoidal, L, K, E);
    const scalar etaMS = etaMeanSq(H, m, T);

    scalar u =
        c*etaCN/h
      - c*(etaCN*etaCN/h/h + etaMS*etaMS/h/h)
     + 1.0/2.0*c*h*(1.0/3.0 - z*z/h/h)*etaX[1];

    scalar w =
       -c*z*(etaX[0]/h*(1.0 - 2.0*etaCN/h) + 1.0/6.0*h*(1.0 - z*z/h/h)*etaX[2]);

    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::waveModels::cnoidal::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    const scalar waveK = mathematical::twoPi/waveLength_;
    const scalar waveKx = waveK*cos(waveAngle_);
    const scalar waveKy = waveK*sin(waveAngle_);

    forAll(level, paddlei)
    {
        const scalar eta =
            this->eta
            (
                waveHeight_,
                m_,
                waveKx,
                waveKy,
                wavePeriod_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                t
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


void Foam::waveModels::cnoidal::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
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
                waveHeight_,
                waterDepthRef_,
                m_,
                waveKx,
                waveKy,
                wavePeriod_,
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

Foam::waveModels::cnoidal::cnoidal
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    regularWaveModel(dict, mesh, patch, false),
    m_(0)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::cnoidal::readDict(const dictionary& overrideDict)
{
    if (regularWaveModel::readDict(overrideDict))
    {
        // Initialise m parameter and wavelength
        initialise
        (
            waveHeight_,
            waterDepthRef_,
            wavePeriod_,
            m_,
            waveLength_
        );

        return true;
    }

    return false;
}


void Foam::waveModels::cnoidal::info(Ostream& os) const
{
    regularWaveModel::info(os);

    os  << "    Cnoidal m parameter : " << m_ << nl;
}


// ************************************************************************* //
