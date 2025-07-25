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

inline Foam::pointVectorField& Foam::pLaplacianMotionSolver::pointMotionU()
{
    return pointMotionU_;
}


inline const Foam::pointVectorField&
Foam::pLaplacianMotionSolver::pointMotionU() const
{
    return pointMotionU_;
}


inline Foam::volVectorField& Foam::pLaplacianMotionSolver::cellMotionU()
{
    return cellMotionU_;
}


inline const Foam::volVectorField&
Foam::pLaplacianMotionSolver::cellMotionU() const
{
    return cellMotionU_;
}


// ************************************************************************* //


