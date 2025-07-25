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

#include "sigmoidalHeaviside.H"
#include "mathematicalConstants.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sigmoidalHeaviside, 1);
addToRunTimeSelectionTable
(
    topOInterpolationFunction,
    sigmoidalHeaviside,
    dictionary
);

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

scalar sigmoidalHeaviside::computeNearBandWidth() const
{
    scalar averageVol(gAverage(mesh_.V().field()));
    const Vector<label>& solutionD = mesh_.solutionD();
    const boundBox& bounds = mesh_.bounds();
    forAll(solutionD, idir)
    {
        if (solutionD[idir] == -1)
        {
            averageVol /= bounds.span()[idir];
            break;
        }
    }

    scalar width = pow(averageVol, scalar(1)/scalar(mesh_.nGeometricD()));
    scalar multMeanRadius =
        dict_.getOrDefaultCompat<scalar>
        (
            "meanRadiusMult", {{"scale", 2306}}, 2
        );
    DebugInfo
        << "Computed near-band width :: " << width
        << " and multiplying with " << multMeanRadius << endl;

    return multMeanRadius*width;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

sigmoidalHeaviside::sigmoidalHeaviside
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    topOInterpolationFunction(mesh, dict),
    dNB_(dict.getOrDefault<scalar>("d", computeNearBandWidth()))
    //dNB_(Function1<scalar>::New("d", dict))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void sigmoidalHeaviside::interpolate
(
    const scalarField& arg,
    scalarField& res
) const
{
    const scalar timeValue = mesh_.time().timeOutputValue();
    const scalar t(timeValue == 0 ? 1. : timeValue);
    const scalar pi = constant::mathematical::pi;
    DebugInfo
        << type() << "::interpolate:: t, dNB " << t << ", " << dNB_ << endl;

    res = 0.5*(scalar(1) + arg/dNB_ + sin(pi*arg/dNB_)/pi);
    res = max(min(scalar(1), res), scalar(0));
}


tmp<scalarField> sigmoidalHeaviside::derivative(const scalarField& arg) const
{
    tmp<scalarField> tderiv = tmp<scalarField>::New(arg.size(), Zero);
    scalarField& deriv = tderiv.ref();
    const scalar timeValue = mesh_.time().timeOutputValue();
    const scalar t(timeValue == 0 ? 1. : timeValue);
    const scalar pi = constant::mathematical::pi;
    scalarField argLimited(max(min(dNB_, arg), -dNB_));
    DebugInfo
        << type() << "::derivative:: t, dNB " << t << ", " << dNB_ << endl;

    deriv = 0.5*(scalar(1) + cos(pi*argLimited/dNB_))/dNB_;

    return tderiv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
