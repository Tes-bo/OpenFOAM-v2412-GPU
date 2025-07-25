/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 IH-Cantabria
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "McCowanWaveModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(McCowan, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        McCowan,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::McCowan::eta
(
    const scalar H,
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar theta,
    const scalar t,
    const scalar X0
) const
{
    vector vec = this->mn(H, h);
    scalar mm = vec[0];
    scalar nn = vec[1];

    scalar C = sqrt(((mag(g_)*h)/mm)*tan(mm));
    scalar ts = 3.5*h/sqrt(H/h);
    scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);

    scalar xin = 0.5*H;
    scalar etas = newtonRapsonF2(xin, H, h, Xa, mm, nn);
    return etas;
}


Foam::vector Foam::waveModels::McCowan::mn
(
    const scalar H,
    const scalar h
) const
{
    // m
    scalar xin = 1;
    scalar m = newtonRapsonF1(xin, H, h);

    // n
    scalar c1 = sin(m + (1.0 + (2.0*H/(3.0*h))));
    scalar n = (2.0/3.0)*sqr(c1);

    return vector(m, n, n);
}


Foam::scalar Foam::waveModels::McCowan::newtonRapsonF1
(
    const scalar x0,
    const scalar H,
    const scalar h
) const
{
    label N = 10000;
    scalar eps = 1.e-5;
    scalar maxval = 10000.0;

    label iter = 1;
    scalar x = x0;
    scalar residual = 0;
    while (iter <= N)
    {
        // f
        scalar a = x + 1.0 + 2.0*H/(3.0*h);
        scalar b = 0.5*x*(1.0 + H/h);
        scalar c = 0.5*x*(1.0 + h/H);
        scalar c1 = sin(a);
        scalar fx = (2.0/3.0)*sqr(c1) - x*H/(h*tan(b));

        residual = mag(fx);

        if (residual < eps)
        {
            return x;
        }
        else if ((iter > 1) && (residual > maxval))
        {
            FatalErrorInFunction
                << "Newton-Raphson iterations diverging: "
                << "iterations = " << iter
                << ", residual = " << residual
                << exit(FatalError);
        }

        // f-prime
        scalar c2 = 1.0/tan(c);
        scalar c3 = 1.0/sin(b);

        scalar fprime = (4.0/3.0)*c1*cos(a) - c2*h/H - b*sqr(c3);

        x -= fx/fprime;
        iter++;
    }

    WarningInFunction
        << "Failed to converge in " << iter << " iterations.  Residual = "
        << residual << nl << endl;

    return x;
}


Foam::scalar Foam::waveModels::McCowan::newtonRapsonF2
(
    const scalar x0,
    const scalar H,
    const scalar h,
    const scalar xa,
    const scalar m,
    const scalar n
) const
{
    label N = 10000;
    scalar eps = 1.e-5;
    scalar maxval = 10000;

    label iter = 1;
    scalar x = x0;
    scalar residual = 0;
    while (iter <= N)
    {
        // f
        scalar a = m*(1.0 + x/h);
        scalar c1 = cos(a);
        scalar c2 = sin(a);

        scalar fx = x - (h*n/m*(c2/(c1 + cosh(m*xa/h))));

        residual = mag(fx);

        if (residual < eps)
        {
            return x;
        }
        else if ((iter > 1) && (residual > maxval))
        {
            FatalErrorInFunction
                << "Newton-Raphson iterations diverging: "
                << "iterations = " << iter
                << ", residual = " << residual
                << exit(FatalError);
        }

        // f-prime
        scalar c3 = cosh(xa*m/h) + c1;
        scalar fprime = 1 - n/c3*(c1 - sqr(c2)/c3);

        x -= fx/fprime;
        iter++;
    }

    WarningInFunction
        << "Failed to converge in " << iter << " iterations.  Residual = "
        << residual << nl << endl;

    return x;
}


Foam::vector Foam::waveModels::McCowan::Uf
(
    const scalar H,
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar theta,
    const scalar t,
    const scalar X0,
    const scalar z
) const
{
    const vector vec = this->mn(H, h);
    const scalar mm = vec[0];
    const scalar nn = vec[1];

    const scalar C = sqrt((mag(g_)*h)/mm*tan(mm));
    const scalar ts = 3.5*h/sqrt(H/h);
    const scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);

    scalar outa = C*nn*(1.0 + cos(mm*z/h)*cosh(mm*Xa/h));
    scalar outb = sqr(cos(mm*z/h) + cosh(mm*Xa/h));

    scalar u = outa/outb;

    outa = C*nn*sin(mm*z/h)*sinh(mm*Xa/h);

    const scalar w = outa/outb;

    const scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


void Foam::waveModels::McCowan::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    forAll(level, paddlei)
    {
        const scalar eta =
            this->eta
            (
                waveHeight_,
                waterDepthRef_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveAngle_,
                t,
                x0_
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::McCowan::McCowan
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    solitaryWaveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::McCowan::readDict(const dictionary& overrideDict)
{
    if (solitaryWaveModel::readDict(overrideDict))
    {
        return true;
    }

    return false;
}


void Foam::waveModels::McCowan::setVelocity
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
                waveHeight_,
                waterDepthRef_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveAngle_,
                t,
                x0_,
                z
            );

            U_[facei] = fraction*Uf*tCoeff;
        }
    }
}


void Foam::waveModels::McCowan::info(Ostream& os) const
{
    solitaryWaveModel::info(os);
}


// ************************************************************************* //
