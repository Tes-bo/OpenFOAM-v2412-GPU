/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "quadratic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quadratic, 0);
    addToRunTimeSelectionTable
    (
        stepUpdate,
        quadratic,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quadratic::quadratic(const dictionary& dict)
:
    stepUpdate(dict),
    minRatio_(coeffsDict().getOrDefault<scalar>("minRatio", 0.1)),
    firstMeritValue_(Zero),
    secondMeritValue_(Zero),
    meritDerivative_(Zero)
{}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::quadratic::updateStep(scalar& step)
{
    Info<< "First merit value,  f(0)  = " << firstMeritValue_ << endl;
    Info<< "Second merit value, f(a0) = " << secondMeritValue_ << endl;
    Info<< "Merit derivative,   df(0) = " << meritDerivative_ << endl;
    Info<< "Previous step,      a0    = " <<  step << endl;
    scalar denom = 1./(step*step);
    scalar coeff1 =
        (secondMeritValue_ - meritDerivative_*step - firstMeritValue_)
      * denom;
    scalar tempStep = - 0.5*meritDerivative_/coeff1;
    if (tempStep < minRatio_*step)
    {
        step = minRatio_*step;
    }
    else
    {
        step = tempStep;
    }
}


void Foam::quadratic::setDeriv(const scalar deriv)
{
    meritDerivative_ = deriv;
}


void Foam::quadratic::setNewMeritValue(const scalar value)
{
    secondMeritValue_ = value;
}


void Foam::quadratic::setOldMeritValue(const scalar value)
{
    firstMeritValue_ = value;
}


// ************************************************************************* //
