/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "wallFunctionCoefficients.H"
#include "dictionary.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Estimate the y+ at the intersection of the two sublayers
static scalar calcYPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11;

    for (int iter = 0; iter < 10; ++iter)
    {
        ypl = log(max(E*ypl, scalar(1)))/kappa;
    }

    return ypl;
}

}  // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallFunctionCoefficients::wallFunctionCoefficients()
:
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{}


Foam::wallFunctionCoefficients::wallFunctionCoefficients
(
    const dictionary& dict
)
:
    Cmu_(dict.getOrDefault<scalar>("Cmu", 0.09)),
    kappa_
    (
        dict.getCheckOrDefault<scalar>("kappa", 0.41, scalarMinMax::ge(SMALL))
    ),
    E_(dict.getCheckOrDefault<scalar>("E", 9.8, scalarMinMax::ge(SMALL))),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallFunctionCoefficients::writeEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<scalar>("Cmu", 0.09, Cmu_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("E", 9.8, E_);
}


// ************************************************************************* //
