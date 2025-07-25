/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    computeSensitivities

Description
    Computes the sensitivities wrt what is defined in the optimisationDict

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "optimisationManager.H"
#include "primalSolver.H"
#include "adjointSolver.H"
#include "incompressibleVars.H"
#include "incompressibleAdjointVars.H"
#include "adjointBoundaryCondition.H"
#include "adjointSolverManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    forAll(adjointSolverManagers, amI)
    {
        PtrList<adjointSolver>& adjointSolvers =
            adjointSolverManagers[amI].adjointSolvers();
        forAll(adjointSolvers, asI)
        {
            adjointSolvers[asI].getObjectiveManager().updateAndWrite();
            adjointSolvers[asI].
                computeObjectiveSensitivities(om.getDesignVariables());
        }
    }

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
