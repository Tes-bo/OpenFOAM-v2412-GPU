/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "SIMPLEControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SIMPLEControl, 0);
    defineRunTimeSelectionTable(SIMPLEControl, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SIMPLEControl::SIMPLEControl
(
    fvMesh& mesh,
    const word& managerType,
    const solver& solver
)
:
    solverControl(solver),
    simpleControl(mesh, "SIMPLE", false),
    managerType_(managerType),
    nIters_(Zero),
    pRefCell_(Zero),
    pRefValue_(Zero)
{
    read();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::SIMPLEControl> Foam::SIMPLEControl::New
(
    fvMesh& mesh,
    const word& managerType,
    const solver& solver
)
{
    auto* ctorPtr = dictionaryConstructorTable(managerType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "control",
            managerType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<SIMPLEControl>(ctorPtr(mesh, managerType, solver));
}


bool Foam::SIMPLEControl::read()
{
    simpleControl::read();
    solverControl::read();
    readIters();

    if (average_ && averageStartIter_ > nIters_)
    {
        WarningInFunction
            << "Average start iteration is larger than nIter in solver "
            << solver_.solverName() << nl
            << tab << "Disabling averaging ..." << nl
            << endl;
        average_ = false;
    }

    return true;
}


void Foam::SIMPLEControl::readIters()
{
    nIters_ = dict().get<label>("nIters");
}


void Foam::SIMPLEControl::checkMeanSolution() const
{
    if (average_ && iter_ < averageStartIter_)
    {
        WarningInFunction
            << "Solver " << solver_.solverName()
            << " converged before averaging started" << nl << tab
            << "Using instantaneous fields ..." << nl
            << endl;
    }
}


bool Foam::SIMPLEControl::converged()
{
    return simpleControl::criteriaSatisfied();
}


bool Foam::SIMPLEControl::performIter()
{
    return true;
}


void Foam::SIMPLEControl::incrementIter()
{
    // Does nothing in base
}


// ************************************************************************* //
