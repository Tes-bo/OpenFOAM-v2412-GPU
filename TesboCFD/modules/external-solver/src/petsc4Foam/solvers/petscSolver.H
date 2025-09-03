/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
    Copyright (C) 2019-2020 Simone Bna
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

Class
    Foam::petscSolver

Group
    grpLduMatrixSolvers

Description
    PETSc solver interface for lduMatrices using a run-time selectable solver
    and preconditioner. PETSc solvers are configured via sub-dictionaries.
    For example,

    \verbatim
    solver petsc;

    petsc
    {
        use_petsc_residual_norm true;

        options
        {
            ksp_type cg;
            ksp_cg_single_reduction true;
        }
    }
    \endverbatim

SourceFiles
    petscSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscSolver_H
#define Foam_module_petscSolver_H

#include "lduMatrix.H"

// PETSc
#include "petscmat.h"
#include "petscksp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class petscSolver Declaration
\*---------------------------------------------------------------------------*/

class petscSolver
:
    public lduMatrix::solver
{
    // Private Data

        //- Coefficients for the petsc solvers
        dictionary petscDict_;

        //- Equation name for solver context
        word eqName_;

        //- Solution name used as prefix
        word prefix_;


    // Private Member Functions

        //- Compute allocation (proc and off-proc)
        //- for the petsc matrix Amat_
        void buildMat
        (
            Mat& Amat,
            List<PetscInt>& lowNonZero,
            PetscInt& maxLowNonZeroPerRow
        ) const;

        //- Update Amat_ matrix by inserting new values
        void updateMat
        (
            Mat& Amat,
            const List<PetscInt>& lowNonZero,
            const PetscInt maxLowNonZeroPerRow
        ) const;

        //- Build the ksp_ krylov solver context
        void buildKsp(Mat& Amat, KSP& ksp) const;

        //- Update the ksp_ krylov solver context
        void updateKsp(KSP& ksp, Mat& Amat, const bool reuse) const;

        //- No copy construct
        petscSolver(const petscSolver&) = delete;

        //- No copy assignment
        void operator=(const petscSolver&) = delete;


public:

    //- Runtime type information
    TypeName("petsc");


    // Constructors

        //- Construct from matrix components and solver controls
        petscSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~petscSolver() = default;


    // Member Functions

        //- Solve with given field and rhs (in solveScalar precision).
        virtual solverPerformance scalarSolve
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt = 0
        ) const;

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
