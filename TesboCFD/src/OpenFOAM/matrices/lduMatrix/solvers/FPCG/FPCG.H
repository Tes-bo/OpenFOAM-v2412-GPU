/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2023 Huawei (Yu Ankun)
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
    Foam::FPCG

Group
    grpLduMatrixSolvers

Description
    A 'faster' preconditioned conjugate gradient solver for
    symmetric lduMatrices using a run-time selectable preconditioner.

    This is termed \em "faster" than the regular PCG since it combines
    global reductions.

SourceFiles
    FPCG.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_FPCG_H
#define Foam_FPCG_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class FPCG Declaration
\*---------------------------------------------------------------------------*/

class FPCG
:
    public lduMatrix::solver
{
    // Private Member Data

        //- Cached preconditioner
        mutable autoPtr<lduMatrix::preconditioner> preconPtr_;


    // Private Member Functions

        //- Blocking version of sum(a*b), sum(mag(b))
        static void gSumMagProd
        (
            FixedList<solveScalar, 2>& globalSum,
            const solveScalarField& a,
            const solveScalarField& b,
            const label comm
        );

        //- No copy construct
        FPCG(const FPCG&) = delete;

        //- No copy assignment
        void operator=(const FPCG&) = delete;


public:

    //- Runtime type information
    TypeName("FPCG");


    // Constructors

        //- Construct from matrix components and solver controls
        FPCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~FPCG() = default;


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance scalarSolve
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt=0
        ) const;

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
