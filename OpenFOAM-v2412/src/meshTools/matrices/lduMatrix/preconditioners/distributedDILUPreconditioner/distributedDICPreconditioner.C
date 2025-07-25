/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "distributedDICPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedDICPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<distributedDICPreconditioner>
        adddistributedDICPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distributedDICPreconditioner::forwardInternalDiag
(
    solveScalarField& rD,
    const label colouri
) const
{
    // Add forward constributions to diagonal

    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    const label nFaces = matrix.upper().size();
    if (cellColourPtr_)
    {
        const auto& cellColour = *cellColourPtr_;
        for (label face=0; face<nFaces; face++)
        {
            const label cell = lPtr[face];
            if (cellColour[cell] == colouri)
            {
                rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[cell];
            }
        }
    }
    else
    {
        for (label face=0; face<nFaces; face++)
        {
            rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[lPtr[face]];
        }
    }
}


void Foam::distributedDICPreconditioner::forwardInternal
(
    solveScalarField& wA,
    const label colouri
) const
{
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    solveScalar* __restrict__ wAPtr = wA.begin();
    const solveScalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    const label nFaces = matrix.upper().size();
    if (cellColourPtr_)
    {
        const auto& cellColour = *cellColourPtr_;
        for (label face=0; face<nFaces; face++)
        {
            const label cell = lPtr[face];
            if (cellColour[cell] == colouri)
            {
                wAPtr[uPtr[face]] -=
                    rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[cell];
            }
        }
    }
    else
    {
        for (label face=0; face<nFaces; face++)
        {
            wAPtr[uPtr[face]] -=
                rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedDICPreconditioner::distributedDICPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& dict
)
:
    distributedDILUPreconditioner(sol, dict)
{}


// ************************************************************************* //
