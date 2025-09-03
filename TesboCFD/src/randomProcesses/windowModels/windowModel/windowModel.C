/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "windowModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(windowModel, 0);
    defineRunTimeSelectionTable(windowModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::windowModel::windowModel(const dictionary& dict, const label nSamples)
:
    scalarField(nSamples),
    overlapPercent_(dict.get<scalar>("overlapPercent")),
    nOverlapSamples_(floor(overlapPercent_/scalar(100)*nSamples)),
    nWindow_(dict.getOrDefault("nWindow", -1))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::windowModel::nSamples() const
{
    return size();
}


Foam::scalar Foam::windowModel::overlapPercent() const
{
    return overlapPercent_;
}


Foam::label Foam::windowModel::nOverlapSamples() const
{
    return nOverlapSamples_;
}


Foam::label Foam::windowModel::nWindow() const
{
    return nWindow_;
}


Foam::label Foam::windowModel::nWindowsTotal(label nSamplesTotal) const
{
    const label nSamples = this->nSamples();

    return floor((nSamplesTotal - nSamples)/(nSamples - nOverlapSamples_)) + 1;
}


Foam::label Foam::windowModel::validate(const label nSamplesTotal)
{
    label nSamples = this->nSamples();

    if (nSamplesTotal < nSamples)
    {
        FatalErrorInFunction
            << "Block size N = " << nSamples
            << " is larger than total number of data points = " << nSamplesTotal
            << exit(FatalError);
    }

    const label nWindowAvailable = nWindowsTotal(nSamplesTotal);

    if (nWindow_ == -1)
    {
        nWindow_ = nWindowAvailable;
    }

    if (nWindow_ > nWindowAvailable)
    {
        FatalErrorInFunction
            << "Number of data points calculated with " << nWindow_
            << " windows greater than the total number of data points"
            << nl
            << "    Block size, N = " << nSamples << nl
            << "    Total number of data points = " << nSamplesTotal << nl
            << "    Maximum number of windows = " << nWindowAvailable << nl
            << "    Requested number of windows = " << nWindow_
            << exit(FatalError);
    }

    const label nRequiredSamples =
        nWindow_*nSamples - (nWindow_ - 1)*nOverlapSamples_;

    Info<< "Windowing:" << nl
        << "    Total samples              : " << nSamplesTotal << nl
        << "    Samples per window         : " << nSamples << nl
        << "    Number of windows          : " << nWindow_ << nl
        << "    Overlap size               : " << nOverlapSamples_ << nl
        << "    Required number of samples : " << nRequiredSamples
        << endl;

    return nRequiredSamples;
}


// ************************************************************************* //
