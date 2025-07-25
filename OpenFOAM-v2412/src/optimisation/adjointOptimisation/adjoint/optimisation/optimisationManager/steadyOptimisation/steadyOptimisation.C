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

#include "steadyOptimisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyOptimisation, 0);
    addToRunTimeSelectionTable
    (
        optimisationManager,
        steadyOptimisation,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyOptimisation::steadyOptimisation(fvMesh& mesh)
:
    optimisationManager(mesh)
{
    initialize();
    // Force the construction of the design variables, if they are not
    // already present
    if (!designVars_)
    {
        designVars_.reset
        (
            designVariables::New
            (
                mesh_,
                subDict("optimisation").subDict("designVariables")
            ).ptr()
        );
    }
    designVars_().setInitialValues();
    dvUpdate_.reset
    (
        new designVariablesUpdate
        (
            mesh,
            subDict("optimisation"),
            adjointSolverManagers_,
            designVars_
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::optimisationManager& Foam::steadyOptimisation::operator++()
{
    time_++;
    if (!end())
    {
        Info<< "\n* * * * * * * * * * * * * * * * *" << endl;
        Info<< "Optimisation cycle " << time_.value() << endl;
        Info<< "* * * * * * * * * * * * * * * * *\n" << endl;
    }
    return *this;
}


Foam::optimisationManager& Foam::steadyOptimisation::operator++(int)
{
    return operator++();
}


bool Foam::steadyOptimisation::checkEndOfLoopAndUpdate()
{
    if (update())
    {
        dvUpdate_->update();
    }
    return end();
}


bool Foam::steadyOptimisation::end()
{
    return time_.end();
}


bool Foam::steadyOptimisation::update()
{
    return (time_.timeIndex() != 1 && !end());
}


// ************************************************************************* //
