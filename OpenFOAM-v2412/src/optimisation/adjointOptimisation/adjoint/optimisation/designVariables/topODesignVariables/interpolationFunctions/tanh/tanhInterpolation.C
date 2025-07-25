/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "tanhInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tanhInterpolation, 0);
    addToRunTimeSelectionTable
    (
        topOInterpolationFunction,
        tanhInterpolation,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tanhInterpolation::tanhInterpolation
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    topOInterpolationFunction(mesh, dict),
    b_(Function1<scalar>::New("b", dict)),
    eta_(dict.getOrDefault<scalar>("eta", scalar(0.5)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tanhInterpolation::interpolate
(
    const scalarField& arg,
    scalarField& res
) const
{
    const scalar time(mesh_.time().timeOutputValue());
    const scalar t(time == 0 ? 1. : time);
    const scalar b(b_->value(t));
    if (debug > 1)
    {
        Info<< type() << "::interpolate:: t, b value " << t << " " << b
            << " eta " << eta_ << endl;
    }

    res = (tanh(b*eta_) + tanh(b*(arg - eta_)))
         /(tanh(b*eta_) + tanh(b*(scalar(1) - eta_)));
}


Foam::tmp<Foam::scalarField> Foam::tanhInterpolation::derivative
(
    const scalarField& arg
) const
{
    tmp<scalarField> tderiv(tmp<scalarField>::New(arg.size(), Zero));
    scalarField& deriv = tderiv.ref();

    const scalar t(mesh_.time().timeOutputValue());
    const scalar b(b_->value(t));
    DebugInfo
        << type() << "::interpolate:: t, b value " << t << " " << b
        << " eta " << eta_ << endl;
    deriv = b*(scalar(1) - sqr(tanh(b*(arg - eta_))))
           /(tanh(b*eta_) + tanh(b*(scalar(1) - eta_)));

    return tderiv;
}


void Foam::tanhInterpolation::setLengthScaleParam(const scalar lengthScale)
{
    eta_ = lengthScale;
}


// ************************************************************************* //
