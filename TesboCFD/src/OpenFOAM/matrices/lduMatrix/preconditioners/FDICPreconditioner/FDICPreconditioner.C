/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "FDICPreconditioner.H"
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FDICPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<FDICPreconditioner>
        addFDICPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FDICPreconditioner::FDICPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD_(sol.matrix().diag().size()),
    rDuUpper_(sol.matrix().upper().size()),
    rDlUpper_(sol.matrix().upper().size())
{
    const scalarField& diag = sol.matrix().diag();
    std::copy(diag.begin(), diag.end(), rD_.begin());

    solveScalar* __restrict__ rDPtr = rD_.begin();
    solveScalar* __restrict__ rDuUpperPtr = rDuUpper_.begin();
    solveScalar* __restrict__ rDlUpperPtr = rDlUpper_.begin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr =
        solver_.matrix().upper().begin();

    const label nCells = rD_.size();
    const label nFaces = solver_.matrix().upper().size();

    for (label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= sqr(upperPtr[face])/rDPtr[lPtr[face]];
    }

    // Generate reciprocal FDIC
    for (label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }

    for (label face=0; face<nFaces; face++)
    {
        rDuUpperPtr[face] = rDPtr[uPtr[face]]*upperPtr[face];
        rDlUpperPtr[face] = rDPtr[lPtr[face]]*upperPtr[face];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FDICPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction
) const
{
    solveScalar* __restrict__ wAPtr = wA.begin();
    const solveScalar* __restrict__ rAPtr = rA.begin();
    const solveScalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().begin();

    const solveScalar* const __restrict__ rDuUpperPtr = rDuUpper_.begin();
    const solveScalar* const __restrict__ rDlUpperPtr = rDlUpper_.begin();

    const label nCells = wA.size();
    const label nFaces = solver_.matrix().upper().size();
    const label nFacesM1 = nFaces - 1;

    for (label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }

    for (label face=0; face<nFaces; face++)
    {
        wAPtr[uPtr[face]] -= rDuUpperPtr[face]*wAPtr[lPtr[face]];
    }

    for (label face=nFacesM1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -= rDlUpperPtr[face]*wAPtr[uPtr[face]];
    }
}


// ************************************************************************* //
