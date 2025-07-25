/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "Hanning.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace windowModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Hanning, 0);
addToRunTimeSelectionTable(windowModel, Hanning, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Hanning::Hanning(const dictionary& dict, const label nSamples)
:
    windowModel(dict, nSamples),
    symmetric_(dict.get<bool>("symmetric")),
    extended_(dict.get<bool>("extended")),
    alpha_(dict.getOrDefault<scalar>("alpha", 0.5))  // Hamming = 0.54
{
    // Extend range if required
    const label offset = extended_ ? 1 : 0;
    const scalar m = nSamples - 1 + 2*offset;

    scalarField t(nSamples);
    forAll(t, i)
    {
        t[i] = i + offset;
    }

    scalarField& wf = *this;
    wf = alpha_ - (1 - alpha_)*cos(constant::mathematical::twoPi*t/m);

    // Reset second half of window if symmetric
    if (symmetric_)
    {
        label midPointI = 0;
        if (nSamples % 2 == 0)
        {
            midPointI = nSamples/2;
        }
        else
        {
            midPointI = (nSamples + 1)/2;
        }

        for (label i = 0; i < midPointI; i++)
        {
            wf[nSamples - i - 1] = wf[i];
        }
    }

    const scalar sumSqr = sum(sqr(wf));

    // Normalisation
    wf *= sqrt(nSamples/sumSqr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Hanning::symmetric() const
{
    return symmetric_;
}


bool Hanning::extended() const
{
    return extended_;
}


Foam::scalar Hanning::alpha() const
{
    return alpha_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace windowModels
} // End namespace Foam

// ************************************************************************* //
