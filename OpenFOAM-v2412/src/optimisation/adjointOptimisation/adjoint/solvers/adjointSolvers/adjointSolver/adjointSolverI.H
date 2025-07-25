/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::adjointSolver::primalSolverName() const
{
    return primalSolverName_;
}


inline const Foam::primalSolver& Foam::adjointSolver::getPrimalSolver() const
{
    return mesh_.lookupObject<primalSolver>(primalSolverName_);
}


inline Foam::primalSolver& Foam::adjointSolver::getPrimalSolver()
{
    return
        const_cast<primalSolver&>
        (
            mesh_.lookupObject<primalSolver>(primalSolverName_)
        );
}


inline const Foam::objectiveManager&
Foam::adjointSolver::getObjectiveManager() const
{
    return objectiveManager_;
}


inline Foam::objectiveManager& Foam::adjointSolver::getObjectiveManager()
{
    return objectiveManager_;
}


inline bool Foam::adjointSolver::isConstraint()
{
    return isConstraint_;
}


inline bool Foam::adjointSolver::isDoubleSidedConstraint()
{
    return isDoubleSidedConstraint_;
}


// ************************************************************************* //
