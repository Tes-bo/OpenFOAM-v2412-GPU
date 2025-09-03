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

#include "lduMatrix.H"
#include "parProfilingSolver.H"
#include "profilingPstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parProfilingSolver, 0);

    typedef lduMatrix::solver baseType;

    addNamedToRunTimeSelectionTable
    (
        baseType,
        parProfilingSolver,
        symMatrix,
        parProfiling
    );

    addNamedToRunTimeSelectionTable
    (
        baseType,
        parProfilingSolver,
        asymMatrix,
        parProfiling
    );
}


// Has been initialised
static bool initialised_(false);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parProfilingSolver::parProfilingSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{
    if (!initialised_)
    {
        initialised_ = true;
        profilingPstream::reset();
        profilingPstream::suspend();
    }

    const word baseSolver(solverControls.get<word>("baseSolver"));

    solvePtr_.reset
    (
        lduMatrix::solver::New
        (
            baseSolver,
            fieldName,
            matrix,
            interfaceBouCoeffs,
            interfaceIntCoeffs,
            interfaces,
            solverControls
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::parProfilingSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    profilingPstream::enable();
    Foam::solverPerformance perf(solvePtr_->solve(psi_s, source, cmpt));
    profilingPstream::suspend();

    return perf;
}


// ************************************************************************* //
