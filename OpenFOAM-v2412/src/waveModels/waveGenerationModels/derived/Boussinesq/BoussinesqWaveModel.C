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

#include "BoussinesqWaveModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Boussinesq, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        Boussinesq,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Boussinesq::eta
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
    scalar C = sqrt(mag(g_)*(H + h));
    scalar ts = 3.5*h/sqrt(H/h);
    scalar aux = sqrt(3.0*H/(4.0*h))/h;
    scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);

    return H*1.0/sqr(cosh(aux*Xa));
}


Foam::vector Foam::waveModels::Boussinesq::Deta
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
    vector deta(Zero);

    scalar C = sqrt(mag(g_)*(H + h));
    scalar ts = 3.5*h/sqrt(H/h);
    scalar a = sqrt(3*H/(4*h))/h;
    scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    scalar expTerm = exp(2*a*Xa);
    scalar b = 8*a*h*expTerm;

    deta[0] =
        b*(1 - expTerm)
       /pow3(1 + expTerm);

    deta[1] =
        2*a*b*(exp(4*a*Xa) - 4*expTerm + 1)
       /pow4(1 + expTerm);

    deta[2] =
       -4*sqr(a)*b*(exp(6*a*Xa) - 11*exp(4*a*Xa) + 11*expTerm - 1)
       /pow5(1 + expTerm);

    return deta;
}


Foam::vector Foam::waveModels::Boussinesq::Uf
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
    scalar C = sqrt(mag(g_)*(H + h));
    scalar eta = this->eta(H, h, x, y, theta, t, X0);
    vector Deta = this->Deta(H, h, x, y, theta, t, X0);

    scalar u =
        C*eta/h
       *(
            1.0
          - eta/(4.0*h)
          + sqr(h)/(3.0*eta)*(1.0 - 3.0/2.0*sqr(z/h))*Deta[1]
        );

    scalar w =
       -C*z/h
       *(
            (1.0 - eta/(2.0*h))*Deta[0]
          + sqr(h)/3.0*(1.0 - 1.0/2.0*sqr(z/h))*Deta[2]
        );

    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::waveModels::Boussinesq::setLevel
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


void Foam::waveModels::Boussinesq::setVelocity
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Boussinesq::Boussinesq
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

bool Foam::waveModels::Boussinesq::readDict(const dictionary& overrideDict)
{
    if (solitaryWaveModel::readDict(overrideDict))
    {
        return true;
    }

    return false;
}


void Foam::waveModels::Boussinesq::info(Ostream& os) const
{
    solitaryWaveModel::info(os);
}


// ************************************************************************* //
