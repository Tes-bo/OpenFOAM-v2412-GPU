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

#include "GrimshawWaveModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Grimshaw, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        Grimshaw,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Grimshaw::alfa
(
    const scalar H,
    const scalar h
) const
{
    scalar eps = H/h;

    return sqrt(0.75*eps)*(1.0 - 0.625*eps + (71.0/128.0)*eps*eps);
}


Foam::scalar Foam::waveModels::Grimshaw::eta
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
    const scalar eps = H/h;
    const scalar eps2 = eps*eps;
    const scalar eps3 = eps*eps2;

    const scalar C =
        sqrt(mag(g_)*h)*sqrt(1.0 + eps - 0.05*eps2 - (3.0/70.0)*eps3);

    const scalar ts = 3.5*h/sqrt(H/h);
    const scalar xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    const scalar alfa = this->alfa(H, h);

    const scalar s = (1.0)/(cosh(alfa*(xa/h)));
    const scalar s2 = s*s;
    const scalar q = tanh(alfa*(xa/h));
    const scalar q2 = q*q;

    return
        h
       *(
            eps*s2
          - 0.75*eps2*s2*q2
          + eps3*(0.625*s2*q2 - 1.2625*s2*s2*q2)
        );
}


Foam::vector Foam::waveModels::Grimshaw::Uf
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
    const scalar eps = H/h;
    const scalar eps2 = eps*eps;
    const scalar eps3 = eps*eps2;

    const scalar C =
        sqrt(mag(g_)*h)*sqrt(1.0 + eps - 0.05*eps2 - (3.0/70.0)*eps3);

    const scalar ts = 3.5*h/sqrt(eps);
    const scalar xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    const scalar alfa = this->alfa(H, h);

    const scalar s = (1.0)/(cosh(alfa*(xa/h)));
    const scalar s2 = s*s;
    const scalar s4 = s2*s2;
    const scalar s6 = s2*s4;

    const scalar zbyh = z/h;
    const scalar zbyh2 = zbyh*zbyh;
    const scalar zbyh4 = zbyh2*zbyh2;

    scalar outa = eps*s2 - eps2*(-0.25*s2 + s4 + zbyh2*(1.5*s2 - 2.25*s4));
    scalar outb = 0.475*s2 + 0.2*s4 - 1.2*s6;
    scalar outc = zbyh2*(-1.5*s2 - 3.75*s4 + 7.5*s6);
    scalar outd = zbyh4*(-0.375*s2 + (45.0/16.0)*s4 - (45.0/16.0)*s6);

    scalar u = sqrt(mag(g_)*h)*(outa - eps3*(outb + outc + outd));

    outa = eps*s2 - eps2*(0.375*s2 + 2*s4 + zbyh2*(0.5*s2 - 1.5*s4));
    outb = (49.0/640.0)*s2 - 0.85*s4 - 3.6*s6;
    outc = zbyh2*((-13.0/16.0)*s2 -(25.0/16.0)*s4 + 7.5*s6);
    outd = zbyh4*(-0.075*s2 -1.125*s4 - (27.0/16.0)*s6);

    const scalar w = sqrt(mag(g_)*h)*(outa - eps3*(outb + outc + outd));

    const scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


void Foam::waveModels::Grimshaw::setLevel
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

Foam::waveModels::Grimshaw::Grimshaw
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

bool Foam::waveModels::Grimshaw::readDict(const dictionary& overrideDict)
{
    if (solitaryWaveModel::readDict(overrideDict))
    {
        return true;
    }

    return false;
}


void Foam::waveModels::Grimshaw::setVelocity
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


void Foam::waveModels::Grimshaw::info(Ostream& os) const
{
    solitaryWaveModel::info(os);
}


// ************************************************************************* //
