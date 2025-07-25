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

#include "solver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solver, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solver::solver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& solverName
)
:
    localIOdictionary
    (
        IOobject
        (
            solverName,
            mesh.time().timeName(),
            fileName("uniform")/fileName("solvers"),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null // avoid type checking since dictionary is read using the
                   // derived type name and type() will result in "solver" here
    ),
    mesh_(mesh),
    managerType_(managerType),
    dict_(dict),
    solverName_(solverName),
    active_(dict.getOrDefault<bool>("active", true)),
    isMaster_(dict.getOrDefault<bool>("isMaster", true)),
    useSolverNameForFields_
        (dict_.getOrDefault<bool>("useSolverNameForFields", false)),
    vars_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solver::readDict(const dictionary& dict)
{
    dict_ = dict;

    // Note: Slightly dangerous to change active_ while the solver is
    // running. At the very least, this should trigger writing before stopping.
    // Additional problems if we have an adjointSolver corresponding to a
    // constraint. To be revisited
    //active_ = dict.getOrDefault<bool>("active", true);

    return true;
}


void Foam::solver::restoreInitValues()
{
    // Does nothing in the base class
}


void Foam::solver::preLoop()
{
    restoreInitValues();
}


void Foam::solver::postLoop()
{
    // Does nothing in the base class
}


void Foam::solver::addTopOFvOptions() const
{
    // Does nothing in the base class
}


// ************************************************************************* //
