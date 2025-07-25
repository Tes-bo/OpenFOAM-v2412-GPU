/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "singleRun.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleRun, 0);
    addToRunTimeSelectionTable(optimisationManager, singleRun, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleRun::singleRun(fvMesh& mesh)
:
    optimisationManager(mesh),
    cycles_(Zero)
{
    initialize();
    if (designVars_)
    {
        designVars_().setInitialValues();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::optimisationManager& Foam::singleRun::operator++()
{
    cycles_++;
    return *this;
}


Foam::optimisationManager& Foam::singleRun::operator++(int)
{
    return operator++();
}


bool Foam::singleRun::checkEndOfLoopAndUpdate()
{
    return end();
}


bool Foam::singleRun::end()
{
    // Force execution of a single loop
    return cycles_ > 1;
}


bool Foam::singleRun::update()
{
    // No update in singleRun cases
    return false;
}


void Foam::singleRun::updateDesignVariables()
{
    // No update in singleRun cases
}


// ************************************************************************* //
