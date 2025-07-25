/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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

Application
    adjointOptimisation

Description
    An automated adjoint-based optimisation loop. Supports multiple types
    of optimisation (shape, topology etc)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "optimisationManager.H"
#include "primalSolver.H"
#include "adjointSolverManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;

    for (om++; !om.end(); om++)
    {
        if (om.update())
        {
            // Update design variables and solve all primal equations
            om.updateDesignVariables();
        }
        else
        {
            // Solve all primal equations
            om.solvePrimalEquations();

            // Clear sensitivities
            om.clearSensitivities();

            // Solve all adjoint equations
            om.solveAdjointEquations();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
