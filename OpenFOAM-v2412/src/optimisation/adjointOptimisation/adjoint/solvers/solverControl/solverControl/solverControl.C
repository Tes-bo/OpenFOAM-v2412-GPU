/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "solverControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solverControl, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::solverControl::read()
{
    // Read basic entries
    printMaxMags_ = solutionDict().getOrDefault<bool>("printMaxMags", false);

    // Manage averaging
    dictionary averagingDict = solutionDict().subOrEmptyDict("averaging");
    averageStartIter_ = averagingDict.getOrDefault<label>("startIter", -1);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solverControl::solverControl(const solver& solver)
:
    solver_(solver),
    printMaxMags_(true),
    iter_(0),
    averageIter_(solver.getOrDefault<label>("averageIter", 0)),
    averageStartIter_(-1),
    // Non run-time modifiable options read in the constructor only
    storeInitValues_
    (
        solverDict().getOrDefault<bool>("storeInitValues", false)
    ),
    average_
    (
        solutionDict().subOrEmptyDict("averaging").
            getOrDefault<bool>("average", false)
    )
{
    read();
}


// ************************************************************************* //
