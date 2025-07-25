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

#include "StokesVWaveModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(StokesV, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        StokesV,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::StokesV::A11
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    return 1.0/s;
}


Foam::scalar Foam::waveModels::StokesV::A13
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);
    return -sqr(c)*(5*sqr(c) + 1)/(8*pow5(s));
}


Foam::scalar Foam::waveModels::StokesV::A15
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
       -(
            1184*pow(c, 10)
          - 1440*pow(c, 8)
          - 1992*pow6(c)
          + 2641*pow4(c)
          - 249*sqr(c) + 18
        )
       /(1536*pow(s, 11));
}


Foam::scalar Foam::waveModels::StokesV::A22
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    return 3/(8*pow4(s));
}


Foam::scalar Foam::waveModels::StokesV::A24
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (192*pow(c, 8) - 424*pow(c, 6) - 312*pow4(c) + 480*sqr(c) - 17)
       /(768*pow(s, 10));
}


Foam::scalar Foam::waveModels::StokesV::A33
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return (13 - 4*sqr(c))/(64*pow(s, 7));
}


Foam::scalar Foam::waveModels::StokesV::A35
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (
            512*pow(c, 12)
          + 4224*pow(c, 10)
          - 6800*pow(c, 8)
          - 12808*pow(c, 6)
          + 16704.0*pow4(c)
          - 3154*sqr(c)
          + 107
        )
       /(4096*pow(s, 13)*(6*sqr(c) - 1));
}


Foam::scalar Foam::waveModels::StokesV::A44
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (80*pow(c, 6) - 816*pow4(c) + 1338*sqr(c) - 197)
       /(1536*pow(s, 10)*(6*sqr(c) - 1));
}


Foam::scalar Foam::waveModels::StokesV::A55
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
       -(
            2880*pow(c, 10)
          - 72480*pow(c, 8)
          + 324000*pow(c, 6)
          - 432000*pow4(c)
          + 163470*sqr(c)
          - 16245
        )
       /(61440*pow(s, 11)*(6*sqr(c) - 1)*(8*pow4(c) - 11*sqr(c) + 3));
}


Foam::scalar Foam::waveModels::StokesV::B22
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return (2*sqr(c) + 1)*c/(4*pow3(s));
}


Foam::scalar Foam::waveModels::StokesV::B24
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (272*pow(c, 8) - 504*pow(c, 6) - 192*pow4(c) + 322*sqr(c) + 21)*c
       /(384*pow(s, 9));
}


Foam::scalar Foam::waveModels::StokesV::B33
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return (8*pow6(c) + 1)*3/(64*pow6(s));
}


Foam::scalar Foam::waveModels::StokesV::B33k
(
    const scalar h,
    const scalar k
) const // d B33 / d k
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    const scalar sk = h*s;
    const scalar ck = h*c;

    return 9.*pow5(c)*ck/(4*pow6(s)) - (9*(8*pow6(c) + 1))/(32*pow(s, 7))*sk;
}


Foam::scalar Foam::waveModels::StokesV::B35
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (
            88128*pow(c, 14)
          - 208224*pow(c, 12)
          + 70848*pow(c, 10)
          + 54000*pow(c, 8)
          - 21816*pow6(c)
          + 6264*pow4(c)
          - 54*sqr(c)
          - 81
        )
       /(12288*pow(s, 12)*(6*sqr(c) - 1));
}


Foam::scalar Foam::waveModels::StokesV::B35k
(
    const scalar h,
    const scalar k
) const // d B35 / d k
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    const scalar sk = h*s;
    const scalar ck = h*c;

    return
        (
            14*88128*pow(c, 13)*ck
          - 12*208224*pow(c, 11)*ck
          + 10*70848*pow(c, 9)*ck
          + 8*54000.0*pow(c, 7)*ck
          - 6*21816*pow5(c)*ck
          + 4*6264*pow3(c)*ck
          - 2*54*c*ck
        )
       /(12288*pow(s, 12)*(6*sqr(c) - 1))
      - (
            88128*pow(c, 14)
          - 208224*pow(c, 12)
          + 70848*pow(c, 10)
          + 54000*pow(c, 8)
          - 21816*pow6(c)
          + 6264*pow4(c)
          - 54*sqr(c)
          - 81
        )*12
       /(12288*pow(s, 13)*(6*sqr(c) - 1))*sk
      - (
            88128*pow(c,14)
          - 208224*pow(c, 12)
          + 70848*pow(c, 10)
          + 54000*pow(c, 8)
          - 21816*pow6(c)
          + 6264*pow4(c)
          - 54*sqr(c)
          - 81
        )*12*c*ck
       /(12288*pow(s, 12)*sqr(6*sqr(c) - 1));
}


Foam::scalar Foam::waveModels::StokesV::B44
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (
            768*pow(c, 10)
          - 448*pow(c, 8)
          - 48*pow6(c)
          + 48*pow4(c)
          + 106*sqr(c)
          - 21
        )*c
       /(384*pow(s, 9)*(6*sqr(c) - 1));
}


Foam::scalar Foam::waveModels::StokesV::B55
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (
            192000*pow(c, 16)
          - 262720*pow(c, 14)
          + 83680*pow(c, 12)
          + 20160*pow(c, 10)
          - 7280*pow(c, 8)
          + 7160*pow(c, 6)
          - 1800*pow(c, 4)
          - 1050*sqr(c)
          + 225
        )
       /(12288*pow(s, 10)*(6*sqr(c) - 1)*(8*pow4(c) - 11*sqr(c) + 3));
}


Foam::scalar Foam::waveModels::StokesV::B55k
(
    const scalar h,
    const scalar k
) const // d B55 / d k
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    const scalar sk = h*s;
    const scalar ck = h*c;

    return
        (
            16*192000*pow(c, 15)*ck
          - 14*262720*pow(c, 13)*ck
          + 12*83680*pow(c, 11)*ck
          + 10*20160*pow(c, 9)*ck
          - 8*7280*pow(c, 7)*ck
          + 6*7160*pow(c, 5)*ck
          - 4*1800*pow(c, 3)*ck
          - 2*1050*pow(c, 1)*ck
        )
       /(12288*pow(s, 10)*(6*sqr(c) - 1)*(8*pow(c, 4) - 11*sqr(c) + 3))
      - (
            192000*pow(c, 16)
          - 262720*pow(c, 14)
          + 83680*pow(c, 12)
          + 20160*pow(c, 10)
          - 7280*pow(c, 8)
          + 7160*pow(c, 6)
          - 1800*pow(c, 4)
          - 1050*pow(c, 2)
          + 225
        )*10.0
       /(12288*pow(s, 11)*(6*sqr(c) - 1)*(8*pow4(c) - 11*sqr(c) + 3))*sk
      - (
            192000*pow(c, 16)
          - 262720*pow(c, 14)
          + 83680*pow(c, 12)
          + 20160*pow(c,10)
          - 7280*pow(c, 8)
          + 7160*pow(c, 6)
          - 1800*pow(c,4)
          - 1050*pow(c, 2)
          + 225
        )*12*c*ck
        /(12288*pow(s, 10)*sqr(6*sqr(c) - 1)*(8*pow4(c) - 11*sqr(c) + 3))
      - (
            192000*pow(c, 16)
          - 262720*pow(c, 14)
          + 83680*pow(c, 12)
          + 20160*pow(c, 10)
          - 7280*pow(c, 8)
          + 7160*pow(c, 6)
          - 1800*pow(c, 4)
          - 1050*pow(c, 2)
          + 225
        )*(32*pow3(c) - 22*c)*ck
        /(12288*pow(s, 10)*(6*sqr(c) - 1)*sqr(8*pow4(c) - 11*sqr(c) + 3));
}


Foam::scalar Foam::waveModels::StokesV::C1
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return (8*pow4(c) - 8*sqr(c) + 9)/(8*pow4(s));
}


Foam::scalar Foam::waveModels::StokesV::C1k
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    const scalar sk = h*s;
    const scalar ck = h*c;

    return
        (4*8*pow3(c)*ck - 2*8*c*ck)/(8*pow4(s))
      - (8*pow4(c) - 8*sqr(c) + 9)*4*sk/(8*pow5(s));
}


Foam::scalar Foam::waveModels::StokesV::C2
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (
            3840*pow(c, 12)
          - 4096*pow(c, 10)
          + 2592*pow(c, 8)
          - 1008*pow(c, 6)
          + 5944*pow(c, 4)
          - 1830*pow(c, 2)
          + 147
        ) // - 2592
       /(512*pow(s, 10)*(6*sqr(c) - 1));
}


Foam::scalar Foam::waveModels::StokesV::C2k
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    const scalar sk = h*s;
    const scalar ck = h*c;

    return
        (
            12*3840*pow(c, 11)*ck
          - 10*4096*pow(c,9)*ck
          + 8*2592*pow(c, 7)*ck
          - 6*1008*pow(c, 5)*ck
          + 4*5944*pow(c, 3)*ck
          - 2*1830*c*ck
        )
       /(512*pow(s, 10)*(6*sqr(c) - 1))
      - (
            3840*pow(c, 12)
          - 4096*pow(c, 10)
          + 2592*pow(c, 8)
          - 1008*pow(c, 6)
          + 5944*pow(c, 4)
          - 1830*pow(c, 2)
          + 147
        )*10*sk
       /(512*pow(s, 11)*(6*sqr(c) - 1))
      - (
            3840*pow(c, 12)
          - 4096*pow(c, 10)
          + 2592*pow(c, 8)
          - 1008*pow(c, 6)
          + 5944*pow(c, 4)
          - 1830*pow(c, 2)
          + 147
        )*12*c*ck
       /(512*pow(s, 10)*sqr(6*sqr(c) - 1));
}


Foam::scalar Foam::waveModels::StokesV::C3
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return -1/(4*s*c);
}


Foam::scalar Foam::waveModels::StokesV::C4
(
    const scalar h,
    const scalar k
) const
{
    const scalar s = sinh(k*h);
    const scalar c = cosh(k*h);

    return
        (12*pow(c, 8) + 36*pow(c, 6) - 162*pow(c, 4) + 141*sqr(c) - 27)
       /(192*c*pow(s, 9));
}


void Foam::waveModels::StokesV::initialise
(
    const scalar H,
    const scalar d,
    const scalar T,
    scalar& kOut,
    scalar& LambdaOut,
    scalar& f1Out,
    scalar& f2Out
) const
{
    scalar f1 = 1;
    scalar f2 = 1;

    const scalar pi = mathematical::pi;
    scalar k = 2.0*pi/(sqrt(mag(g_)*d)*T);
    scalar lambda = H/2.0*k;

    label n = 0;

    static const scalar tolerance = 1e-12;
    static const label iterMax = 10000;

    while ((mag(f1) > tolerance || mag(f2) > tolerance) && (n < iterMax))
    {
        const scalar b33 = B33(d, k);
        const scalar b35 = B35(d, k);
        const scalar b55 = B55(d, k);
        const scalar c1 = C1(d, k);
        const scalar c2 = C2(d, k);

        const scalar b33k = B33k(d, k);
        const scalar b35k = B35k(d, k);
        const scalar b55k = B55k(d, k);
        const scalar c1k = C1k(d, k);
        const scalar c2k = C2k(d, k);

        const scalar l2 = sqr(lambda);
        const scalar l3 = l2*lambda;
        const scalar l4 = l3*lambda;
        const scalar l5 = l4*lambda;

        const scalar Bmat11 =
            2*pi/(sqr(k)*d)*(lambda + l3*b33 + l5*(b35 + b55))
          - 2*pi/(k*d)*(l3*b33k + l5*(b35k + b55k));

        const scalar Bmat12 =
          - 2*pi/(k*d)*(1 + 3*l2*b33 + 5*l4*(b35 + b55));

        const scalar Bmat21 =
          - d/(2*pi)*tanh(k*d)*(1 + l2*c1 + l4*c2)
          - k*d/(2*pi)*(1 - sqr(tanh(k*d)))*d*(1 + l2*c1 + l4*c2)
          - k*d/(2*pi)*tanh(k*d)*(l2*c1k + l4*c2k);

        const scalar Bmat22 = - k*d/(2.0*pi)*tanh(k*d)*(2*lambda*c1 + 4*l3*c2);

        f1 = pi*H/d - 2*pi/(k*d)*(lambda + l3*b33 + l5*(b35 + b55));

        f2 =
        (
            (2*pi*d)/(mag(g_)*sqr(T))
          - k*d/(2*pi)*tanh(k*d)*(1 + l2*c1 + l4*c2)
        );

        const scalar lambdaPr =
            (f1*Bmat21 - f2*Bmat11)/(Bmat11*Bmat22 - Bmat12*Bmat21);
        const scalar kPr =
            (f2*Bmat12 - f1*Bmat22)/(Bmat11*Bmat22 - Bmat12*Bmat21);

        lambda += lambdaPr;
        k += kPr;

        n++;
    }

    kOut = k;
    LambdaOut = lambda;

    f1Out = mag(f1);
    f2Out = mag(f2);
}


Foam::scalar Foam::waveModels::StokesV::eta
(
    const scalar h,
    const scalar kx,
    const scalar ky,
    const scalar lambda,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar t,
    const scalar phase
) const
{
    const scalar k = sqrt(kx*kx + ky*ky);

    const scalar b22 = B22(h, k);
    const scalar b24 = B24(h, k);
    const scalar b33 = B33(h, k);
    const scalar b35 = B35(h, k);
    const scalar b44 = B44(h, k);
    const scalar b55 = B55(h, k);

    const scalar l2 = sqr(lambda);
    const scalar l3 = l2*lambda;
    const scalar l4 = l3*lambda;
    const scalar l5 = l4*lambda;

    const scalar amp1 = lambda/k;
    const scalar amp2 = (b22*l2 + b24*l4)/k;
    const scalar amp3 = (b33*l3 + b35*l5)/k;
    const scalar amp4 = b44*l4/k;
    const scalar amp5 = b55*l5/k;

    const scalar theta = kx*x + ky*y - 2.0*mathematical::pi/T*t + phase;

    return
    (
        amp1*cos(theta)
      + amp2*cos(2*theta)
      + amp3*cos(3*theta)
      + amp4*cos(4*theta)
      + amp5*cos(5*theta)
    );
}


Foam::vector Foam::waveModels::StokesV::Uf
(
    const scalar d,
    const scalar kx,
    const scalar ky,
    const scalar lambda,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar t,
    const scalar phase,
    const scalar z
) const
{
    const scalar k = sqrt(kx*kx + ky*ky);

    const scalar a11 = A11(d, k);
    const scalar a13 = A13(d, k);
    const scalar a15 = A15(d, k);
    const scalar a22 = A22(d, k);
    const scalar a24 = A24(d, k);
    const scalar a33 = A33(d, k);
    const scalar a35 = A35(d, k);
    const scalar a44 = A44(d, k);
    const scalar a55 = A55(d, k);

    const scalar pi = mathematical::pi;
    const scalar l2 = sqr(lambda);
    const scalar l3 = l2*lambda;
    const scalar l4 = l3*lambda;
    const scalar l5 = l4*lambda;

    const scalar a1u = 2*pi/T/k*(lambda*a11 + l3*a13 + l5*a15);
    const scalar a2u = 2*2*pi/T/k*(l2*a22 + l4*a24);
    const scalar a3u = 3*2*pi/T/k*(l3*a33 + l5*a35);
    const scalar a4u = 4*2*pi/T/k*(l4*a44);
    const scalar a5u = 5*2*pi/T/k*(l5*a55);

    const scalar theta = kx*x + ky*y - 2*pi/T*t + phase;

    scalar u =
        a1u*cosh(k*z)*cos(theta)
      + a2u*cosh(2*k*z)*cos(2*theta)
      + a3u*cosh(3*k*z)*cos(3*theta)
      + a4u*cosh(4*k*z)*cos(4*theta)
      + a5u*cosh(5*k*z)*cos(5*theta);

    scalar w =
        a1u*sinh(k*z)*sin(theta)
      + a2u*sinh(2*k*z)*sin(2*theta)
      + a3u*sinh(3*k*z)*sin(3*theta)
      + a4u*sinh(4*k*z)*sin(4*theta)
      + a5u*sinh(5*k*z)*sin(5*theta);

    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::waveModels::StokesV::setLevel
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
                waterDepthRef_,
                waveKx,
                waveKy,
                lambda_,
                wavePeriod_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                t,
                wavePhase_
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


void Foam::waveModels::StokesV::setVelocity
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
                waterDepthRef_,
                waveKx,
                waveKy,
                lambda_,
                wavePeriod_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                t,
                wavePhase_,
                z
            );

            U_[facei] = fraction*Uf*tCoeff;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::StokesV::StokesV
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    StokesI(dict, mesh, patch, false),
    lambda_(0)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::StokesV::readDict(const dictionary& overrideDict)
{
    if (StokesI::readDict(overrideDict))
    {
        scalar f1;
        scalar f2;
        scalar waveK;

        initialise
        (
            waveHeight_,
            waterDepthRef_,
            wavePeriod_,
            waveK,
            lambda_,
            f1,
            f2
        );

        if (f1 > 0.001 || f2 > 0.001)
        {
            FatalErrorInFunction
                << "No convergence for Stokes V wave theory" << nl
                << "    f1: " << f1 << nl
                << "    f2: " << f2 << nl
                << exit(FatalError);
        }

        return true;
    }

    return false;
}


void Foam::waveModels::StokesV::info(Ostream& os) const
{
    StokesI::info(os);

    os  << "    Lambda : " << lambda_ << nl
        << "    Wave type : " << waveType() << nl;
}


// ************************************************************************* //
