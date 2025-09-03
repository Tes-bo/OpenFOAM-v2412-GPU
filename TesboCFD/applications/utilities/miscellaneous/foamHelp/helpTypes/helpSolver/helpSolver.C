/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "helpSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace helpTypes
    {
        defineTypeNameAndDebug(helpSolver, 0);
        addNamedToRunTimeSelectionTable
        (
            helpType,
            helpSolver,
            dictionary,
            solver
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helpTypes::helpSolver::helpSolver()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helpTypes::helpSolver::~helpSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpTypes::helpSolver::init()
{
    helpType::init();

    argList::addBoolOption
    (
        "read",
        "Read solver type from the system/controlDict"
    );
}


void Foam::helpTypes::helpSolver::execute
(
    const argList& args,
    const fvMesh& mesh
)
{
    word solver;

    if (args.readIfPresent("browse", solver))
    {
        displayDoc(solver, ".*solvers/.*Foam/", true, "C");
    }
    else if (args.found("read"))
    {
        mesh.time().controlDict().readEntry("application", solver);
        displayDoc(solver, ".*solvers/.*Foam/", true, "C");
    }
    else
    {
        displayDocOptions(".*solvers/.*Foam/", true, "C");
    }
}


// ************************************************************************* //
