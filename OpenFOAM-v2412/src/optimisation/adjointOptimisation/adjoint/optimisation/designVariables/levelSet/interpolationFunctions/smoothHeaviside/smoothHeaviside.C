/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "smoothHeaviside.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

defineTypeNameAndDebug(smoothHeaviside, 0);
addToRunTimeSelectionTable
(
    topOInterpolationFunction,
    smoothHeaviside,
    dictionary
);

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

smoothHeaviside::smoothHeaviside
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    topOInterpolationFunction(mesh, dict),
    b_(Function1<scalar>::New("b", dict))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void smoothHeaviside::interpolate
(
    const scalarField& arg,
    scalarField& res
) const
{
    const scalar timeValue = mesh_.time().timeOutputValue();
    const scalar t(timeValue == 0 ? 1. : timeValue);
    const scalar b(b_->value(t));
    res = 0.5*(scalar(1) + tanh(b*arg));
}


tmp<scalarField> smoothHeaviside::derivative(const scalarField& arg) const
{
    tmp<scalarField> tderiv = tmp<scalarField>::New(arg.size(), Zero);
    scalarField& deriv = tderiv.ref();
    const scalar timeValue = mesh_.time().timeOutputValue();
    const scalar t(timeValue == 0 ? 1. : timeValue);
    const scalar b(b_->value(t));

    deriv = 0.5*b*(scalar(1) - sqr(tanh(b*arg)));

    return tderiv;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
