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

#include "ArmijoConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ArmijoConditions, 1);
    addToRunTimeSelectionTable
    (
        lineSearch,
        ArmijoConditions,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ArmijoConditions::ArmijoConditions
(
    const dictionary& dict,
    const Time& time,
    updateMethod& UpdateMethod
)
:
    lineSearch(dict, time, UpdateMethod),
    c1_(coeffsDict().getOrDefault<scalar>("c1", 1.e-4))
{}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool Foam::ArmijoConditions::converged()
{
    DebugInfo
        << "New merit function value " << newMeritValue_ << endl;
    DebugInfo
        << "Old merit function value " << oldMeritValue_ << endl;
    DebugInfo
        << "c1, step, directionalDeriv "
        << c1_ << " " << step_ << " " <<directionalDeriv_
        << endl;
    DebugInfo<< "Extrapolated merit function value "
        << oldMeritValue_ + c1_*step_*directionalDeriv_ << endl;
    return newMeritValue_ < oldMeritValue_ + c1_*step_*directionalDeriv_;
}


void Foam::ArmijoConditions::updateStep()
{
    stepUpdate_->updateStep(step_);
    Info<< "Using step " << step_ << endl;
}


// ************************************************************************* //
