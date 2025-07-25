/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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

const Foam::fvMesh& Foam::solver::mesh() const
{
    return mesh_;
}


const Foam::word& Foam::solver::solverName() const
{
    return solverName_;
}


bool Foam::solver::useSolverNameForFields() const
{
    return useSolverNameForFields_;
}


Foam::word Foam::solver::extendedVariableName(const word& varName) const
{
    return (useSolverNameForFields_ ? word(varName + solverName_) : varName);
}


bool Foam::solver::active()
{
    return active_;
}


const Foam::dictionary& Foam::solver::dict() const
{
    return dict_;
}


const Foam::variablesSet& Foam::solver::getVariablesSet() const
{
    return vars_();
}


Foam::variablesSet& Foam::solver::getVariablesSet()
{
    return vars_();
}


const Foam::word& Foam::solver::managerType() const
{
    return managerType_;
}


bool Foam::solver::isMaster() const
{
    return isMaster_;
}


// ************************************************************************* //
